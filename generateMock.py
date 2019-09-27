import os
from os import path
import gzip
import argparse
import logging
from collections import Counter
import pandas as pd
import numpy as np

def processRunID(sample):
  dirpath = "/k11e/pvdisk/bigbase/kbdata/sampledata"
  sample_path = path.join(dirpath, sample[0:1],sample[0:2], sample[0:3], sample[0:4])
  R1 = path.join(sample_path, sample+'.R1.fastq.gz')
  R2 = path.join(sample_path, sample+'.R2.fastq.gz')
  return {"read1":R1, "read2":R2}


class GenerateMockData:
  '''info is instance of pd.DataFrame, two fields is required,run_id,relAbundance,read1,read2;
  data is the reads number to generate;
  '''
  def __init__(self, info, outdir='./Mock',r1="Mock.R1.fastq.gz", r2="Mock.R2.fastq.gz", data=3333333):
    #info['readsNumber'] = (data*info.relAbundance/100).astype(int)
    choice = np.random.choice(info.index, data, p=info.relAbundance/info.relAbundance.sum())
    info['readsNumber'] = pd.Series(Counter(choice))
    info.readsNumber.fillna(0, downcast='infer', inplace=True)
    self.info = info
    outdir = path.abspath(outdir)
    self.r1 = path.join(outdir, r1)
    self.r2 = path.join(outdir, r2)
  def __enter__(self):
    self.open()
    return self
  def __exit__(self, exc_type, exc_val, exc_tb):
    self.close()
  def open(self):
    self.r1 = gzip.open(self.r1, 'wb')
    self.r2 = gzip.open(self.r2, 'wb')
  def close(self):
    self.r1.close()
    self.r2.close()
  def __iter__(self):
    self.index = 0
    return self
  def __next__(self):
    if self.index < self.sampleN:
      self.processOne()
      self.r1.flush()
      self.r2.flush()
      self.index += 1
    else:
      self.close()
      raise StopIteration
  @property
  def status(self):
    return self.info.iloc[self.index,:]
  @property
  def sampleN(self):
    return self.info.shape[0]
  def processOne(self):
    logger.info("start " + self.status.run_id)
    logger.info(self.status.run_id + "\tAbundance readsNumber: " + str(self.status.readsNumber))
    with gzip.open(self.status.read1, 'rb') as f:
      fq = bytearray(f.read()).strip().split(b'\n')
      logger.debug("bytearray is complete")
    with gzip.open(self.status.read2, 'rb') as f:
      fq2 = bytearray(f.read()).strip().split(b'\n')
    subfix = ('@' + self.status.run_id+'_not_').encode()
    is_replace = True if self.status.readsNumber > len(fq)//4 else False
    randomChr = list(map(lambda x: chr(x), np.arange(65,91)))
    #for i in fq[0::4]: i[0:1] = subfix
    if self.status.readsNumber > 0:
      #selected = np.random.randint(0, len(fq)//4, self.status.readsNumber)
      selected = np.random.choice(np.arange(len(fq)//4), self.status.readsNumber, replace=is_replace)
    else:
      return
    if not is_replace:
      for i in selected:
        fq[4*i][0:1] = subfix
        fq2[4*i][0:1] = subfix
    else:
      for i in selected:
        substitute = subfix.replace(b'not', (''.join(np.random.choice(randomChr, 3))).encode())
        fq[4*i][0:1] = substitute
        fq2[4*i][0:1] = substitute
    logger.debug("run_id add complete")
    list(map(lambda x: self.r1.write(b'\n'.join(fq[(x*4):(x*4+4)])+b'\n'), selected))
    list(map(lambda x: self.r2.write(b'\n'.join(fq2[(x*4):(x*4+4)])+b'\n'), selected))
    logger.info("complete: " + self.status.run_id)
  #def generateReads(self):
    
def parseArgs():
  parser = argparse.ArgumentParser(description='Generate Mock fastq from fastq')
  parser.add_argument("--info", "-i", required=True, help="infomation for generate reads")
  parser.add_argument("--outdir", "-o", required=False, default='.', help="output direction") 
  parser.add_argument("--nsample", "-n", required=False, default=1, type=int, help="the number samples to generate")
  parser.add_argument("--fname", "-f", required=False, default="Mock", help='fastq file name')
  parser.add_argument("--readsNumber", "-r", required=False, default=3333333, type=int, help="reads number to generate")
  parser.add_argument('--log', '-l', required=False, default='INFO', choices=['WARNING', 'INFO', 'DEBUG'], help='log level, WARNING, INFO, DEBUG')
  #parser.add_argument("--mockName", "-m", required=False, default=None, help="one line for one mock")
  args = parser.parse_args() 
  return args

if __name__=='__main__':
  args = parseArgs()
  info = pd.read_csv(args.info, sep='\t', header=0)
  info.rename(columns={'DNA_quality_neq_rel':'relAbundance'}, inplace=True)
  # 'relAbundance', 'read1', 'read2' are required
  keepColumns = ['species', 'genus', 'run_id', 'relAbundance', 'read1', 'read2']
  read = pd.DataFrame(list(info.run_id.apply(processRunID)))
  info["read1"] = read.loc[:, "read1"]
  info["read2"] = read.loc[:, "read2"]
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
    mock = GenerateMockData(info=info.loc[:, keepColumns], outdir=outdir,
                            r1=mockName+".R1.fastq.gz", r2=mockName+".R2.fastq.gz",
                            data=data)
    mock.open()
    iter(mock)
    list(mock)
    mock.close()




