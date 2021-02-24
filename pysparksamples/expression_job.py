# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
### Glue job to fetch and aggregate gene expression data from TCGA
## Command line arguments:
## --project - the TCGA project (e.g TCGA-BRCA for breast cancer)
## --output_bucket - the output bucket where results are to be written to
## --workflow_type - the Workflow type defines the output (default: HTSeq - FPKM-UQ)

import sys
import boto3
import requests 
import json
import pandas as pd
import gzip

from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from pyspark.context import SparkContext
from awsglue.context import GlueContext
from awsglue.job import Job
from pyspark.sql.types import Row, StringType, IntegerType, ArrayType, StructType, DoubleType, BooleanType, DateType
from  pyspark.sql.functions import input_file_name, concat, col
from pyspark.sql.functions import first, last
from pyspark.sql.types import IntegerType
from pyspark.sql.functions import udf, struct, explode

sc = SparkContext.getOrCreate()
glueContext = GlueContext(sc)
spark = glueContext.spark_session
job = Job(glueContext)


## Get the argumentlist

args=getResolvedOptions(sys.argv,
                        ['JOB_NAME',
                        'project',
                        'workflow_type',
                        'output_bucket'])

## The GDC endpoint for files and the NCI endpoint to query for the S3 URL

files_endpt = 'https://api.gdc.cancer.gov/files'
indexd_endpt = 'https://nci-crdc.datacommons.io/index/index/'

s3_tcga_bucket = 'tcga-2-open'
s3 = boto3.resource('s3')
output_bucket= args['output_bucket']

project_id = args['project']


## method to query the NCI endpoint for the S3 path
## Inputs to this method are the UUID and submitter ID from the GDC endpoint query


    
    
## Fields to be returned as a comma separated list

fields = [
      "file_name"
    , "cases.primary_site"
    , "cases.case_id"
    , "cases.project.project_id"
    , "cases.submitter_id"
    , "cases.samples.submitter_id"
    , "cases.samples.sample_id"

]


size = 5000
fields = ','.join(fields)

#data_type = args['data_type']
workflow_type = args['workflow_type']

data_type = "Gene Expression Quantification"
#workflow_type = "HTSeq - FPKM-UQ"

workflow = workflow_type.replace(" ","")

filters = {
    "op":"and",
    "content":[
        {"op": "in",
        "content":{
            "field": "cases.project.project_id",
            "value": [project_id]
            }
        },
        {"op": "in",
        "content":{
            "field": "files.data_type",
            "value": [data_type]
            }
        },
        {"op": "in",
        "content":{
            "field": "files.analysis.workflow_type",
            "value": [workflow_type]
            }
        }
    ]
}
# With a GET request, the filters parameter needs to be converted
# from a dictionary to JSON-formatted string

params = {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "JSON",
    "size": size
    }

## query the files endpoint and get back JSON response

query_response = requests.get(files_endpt, params = params)

json_response = json.loads(query_response.content.decode("utf-8"))["data"]["hits"]

print(len(json_response))
## Parallel read of JSON object

df = spark.read.json(sc.parallelize([json_response]))
df2 = df.repartition(8)
uf = df2.select("id",explode(df2.cases.samples)).toDF("id","samples").select('id','samples.submitter_id','samples.sample_id')

def get_data(uuid, sample_submitter_id):
    query_response = requests.get(indexd_endpt + "/" + uuid)
    urls_response = json.loads(query_response.content.decode("utf-8"))["urls"]
    url = [x for x in urls_response if x.startswith("s3://")]
    if len(url) != 1:
        print("Something weird with UUID " + uuid + "returned " + str(url))
    url = url[0]
    return url
    
urldf = udf(get_data)

## construct list of S3 input paths

inputpath=uf.withColumn('Result', urldf('id', 'submitter_id'))
inputlist = list(inputpath.select('Result').toPandas()['Result'])


## This section parallelizes reading of all the files and stores data in a single data frame

schema = schema = StructType() \
    .add("GeneID", StringType(), True) \
    .add("Expression",DoubleType(), True)

data = spark.read.option("sep","\t").csv(inputlist,header=False,mode="DROPMALFORMED",schema=schema)

## Add a column with file path, which adds the s3 file path from which row is extracted
data = data.withColumn("fullpath", input_file_name())

## Add a column which is a substring of full s3 path and gives filename so that we can match and join with Json data
data = data.withColumn("file", data.fullpath.substr(55,100))
data = data.withColumn("EnsemblGene",data.GeneID.substr(1,15))

## Join the data with the frame that includes submitter ID

joinedf = data.join(inputpath,data["fullpath"]==inputpath["Result"])

## Select only relevant columns

allgene = joinedf.select("EnsemblGene" ,joinedf["submitter_id"].getItem(0), "Expression").withColumnRenamed("submitter_id[0]","submitter_id")

## Pivot the dataframe to form the expression matrix

pivotframe = allgene.groupBy('EnsemblGene').pivot('submitter_id').agg(first("Expression"))

## Write the data frame to the output bucket

pivotframe.write.mode("overwrite").parquet("s3://"+ output_bucket + "/tcga-expression/" + project_id + "/" + workflow)

job.commit()
