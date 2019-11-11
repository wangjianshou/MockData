import sys
import time
import re
import ssl
import requests
from bs4 import BeautifulSoup
ssl._create_default_https_context = ssl._create_unverified_context

class DNAlength:
  def __init__(self, micro):
    self.micro = micro
  def __iter__(self):
    self.index = 0
    return self
  def MicroLength(self, term):
    url = "https://www.ncbi.nlm.nih.gov/genome/?term=" + term
    page = requests.get(url).text
    soup = BeautifulSoup(page, "html5lib")
    title = soup.find("span", "GenomeTitle")
    speciesName = hasattr(title, "text") and title.text or "null" 
    try:
      #GC = soup.find("table","summary").findAll("tr")[4].text.split(': ')[-1]
      #length = soup.find("table","summary").findAll("tr")[2].text.split(': ')[-1]
      length = soup.find('td', text=re.compile(".*total length.*")).text.split(": ")[-1]
      GC = soup.find('td', text=re.compile(".*GC%.*")).text.split(": ")[-1]
    except AttributeError:
      GC = "null"
      length = "null"
    return (term, speciesName, length, GC)
  def __next__(self):
    if self.index >= len(self.micro):
      raise StopIteration
    else:
      r = self.MicroLength(self.micro[self.index])
      self.index += 1
      print(self.index)
      return r

if __name__=="__main__":
  with open(sys.argv[1], "r") as f:
    micro = [i.strip() for i in f]
  obj = DNAlength(micro)
  iter(obj)
  with open("genomeinfo.txt", "w") as f:
    f.write("\t".join(("searchName","ncbiName","Length(Mb)","GC%\n")))
    for i in obj:
      #ssl._create_default_https_context = ssl._create_unverified_context
      time.sleep(10)
      f.write('\t'.join(i)+'\n')
      f.flush()









  
