#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTQ_DIR="$PROJECT_DIR/data/fastq"

mkdir -p "$FASTQ_DIR"

cd "$FASTQ_DIR"

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/022/SRR10416722/SRR10416722_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/021/SRR10416721/SRR10416721_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/016/SRR10416716/SRR10416716_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/019/SRR10416719/SRR10416719_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/012/SRR10416712/SRR10416712_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/015/SRR10416715/SRR10416715_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/014/SRR10416714/SRR10416714_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/013/SRR10416713/SRR10416713_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/020/SRR10416720/SRR10416720_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/021/SRR10416721/SRR10416721_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/018/SRR10416718/SRR10416718_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/022/SRR10416722/SRR10416722_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/011/SRR10416711/SRR10416711_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/016/SRR10416716/SRR10416716_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/017/SRR10416717/SRR10416717_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/013/SRR10416713/SRR10416713_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/020/SRR10416720/SRR10416720_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/014/SRR10416714/SRR10416714_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/018/SRR10416718/SRR10416718_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/017/SRR10416717/SRR10416717_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/011/SRR10416711/SRR10416711_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/012/SRR10416712/SRR10416712_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/015/SRR10416715/SRR10416715_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/019/SRR10416719/SRR10416719_1.fastq.gz
