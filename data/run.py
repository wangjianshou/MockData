import sys
import os
import gzip
import json
import requests
import pandas as pd
from os import path
import argparse
import logging
from collections import Counter
import numpy as np
import pandas as pd

sys.path.append('/k11e/pvdisk/fastbase/Users/wangjianshou/git/MockData')
from src.generateMock import GenerateMockData,parseArgs


args = parseArgs()
info = pd.read_csv(args.info, sep='\t', header=0)
loggerLevel = 10 if args.log=='WARNING' else 20 if args.log=='INFO' else 30
logger = logging.getLogger(__name__)
logging.basicConfig(level=loggerLevel,format = '%(asctime)s - %(levelname)s - %(message)s',
                    filename=path.join(args.outdir, 'log.GenerateMockReads'), filemode='w')

m = args.readsNumber
s = 0.05 * m
for i in range(args.nsample):
  while True:
    data = int(np.random.normal(loc=m, scale=s))
    if data > m-s and data < m+s: break
  mockName = args.fname + '_' + str(i)
  logger.info(mockName+'_readsNumber:'+str(data))
  outdir = path.join(path.abspath(args.outdir), mockName)
  os.path.isdir(outdir) or os.makedirs(outdir)
  mock = GenerateMockData(info=info, logger=logger, outdir=outdir,
                          r1=mockName+".R1.fastq.gz", r2=mockName+".R2.fastq.gz",
                          data=data)
  mock.open()
  iter(mock)
  list(mock)
  mock.close()
