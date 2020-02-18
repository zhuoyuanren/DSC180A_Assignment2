import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import shutil
import shlex
import sys
import subprocess as sp

def filter_and_recode(input_file, output_dir = "data/interim/", file_name="chr22"):
    if not os.path.exists(output_dir):
        cmd = shlex.split("mkdir -p "+output_dir)
        sp.call(cmd)
    cmd = shlex.split("plink2  --vcf "+input_file+"   --make-bed   --snps-only   --maf 0.05   --geno 0.1   --mind 0.05   --recode  --out "+ output_dir+'/'+file_name)
    sp.call(cmd)

def pca(data_dir="data/interim/chr22", remove_file=None):
    if remove_file:
        cmd = shlex.split("plink2 --bfile "+data_dir+" --pca  --remove  "+remove_file+"  --out "+data_dir)
        sp.call(cmd)
    else:
        cmd = shlex.split("plink2 --bfile "+data_dir+" --pca  --out "+data_dir)
        sp.call(cmd)
        
def visualize(data_dir="data/interim/chr22", outlier_file=None):
    eigen = pd.read_table(data_dir+'.eigenvec', header=None, sep=' ')
    eigen= eigen[range(1,12)].set_index([1])    
    ax = eigen.plot(x=2,y=3,kind='scatter')
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    plt.savefig('data/pca_plot.png')
    if outlier_file:
        std = eigen.std().mean()
        outliers = (eigen.abs()<2*std).all(axis=1)
        outlier_lst=list(outliers[outliers==False].index)
        with open(outlier_file,'w') as filehandle:
            for listitem in outlier_lst:
                filehandle.write('%s %s\n' % (listitem, listitem))

parser = argparse.ArgumentParser(description='PCA and visualization with Plink2')
parser.add_argument('data', metavar='data_dir', type=str, nargs=1, help='dir to the .vcf.gz file')
parser.add_argument('process', type=str, nargs=1, help='the process to deal with')
args = parser.parse_args()

print(args.process[0])

if args.process[0]=="filter":
    filter_and_recode(args.data[0])
elif args.process[0]=="pca":
    pca()
elif args.process[0]=="visualize":
    visualize()
elif args.process[0]=="outlier":
    visualize(outlier_file="data/outlier.txt")
elif args.process[0]=="pca_remove":
    pca(remove_file = "data/outlier.txt")