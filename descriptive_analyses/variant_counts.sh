# 2023-06-16
# Gwenna Breton
# Goal: count known and unreported variants in our newly generated data using an updated version of dbsnp.

# Preliminary step: retrieve and prepare dbsnp file
ssh gwennabr@rackham.uppmax.uu.se
screen -S dbsnp
cd /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156
interactive -A p2018003 -p core -n 2 -t 2:00:00 -M snowy
#To list the files:
curl ftp://ftp.ncbi.nih.gov/snp/redesign/archive/b156/VCF/
#To download the files:
wget ftp://ftp.ncbi.nih.gov/snp/redesign/archive/b156/VCF/GCF_000001405.40.gz.md5 #And all the other files in the folder

#The relevant file for hg38 is GCF_000001405.40.gz. However the contigs are named like NC_000001.11 for chr1, etc. I need to change the contig names to what I have in my VCF files.
# I decided to keep only 1-22XYMT (I could skip Y and MT actually) and rename these to what I have in my files.
cd /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/modified_version
zgrep -v "##" ../VCF/GCF_000001405.40.gz | cut -f1 | sort | uniq > list_of_contigs
grep NC list_of_contigs > list_of_contigs_conversion_1-22XYMT
#Then I made a list that listed the "new" names and pasted the two lists together.

#select 1-22XYMT and rename the contigs (trying to do all in one go)

#The following has been written to /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/modified_version/select_rename_dbsnp.sh
sbatch /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/modified_version/select_rename_dbsnp.sh

#!/bin/bash -l
 
#SBATCH -A p2018003
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J select_rename_dbsnp
#SBATCH -o select_rename_dbsnp.output
#SBATCH -e select_rename_dbsnp.output

module load bioinfo-tools
module load bcftools/1.17
invcf=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/VCF/GCF_000001405.40.gz
outvcf=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/modified_version/dbsnp156_1-22XYMT.vcf.gz
conversion=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/modified_version/list_of_contigs_conversion_1-22XYMT # "old_name new_name\n" pairs separated by whitespaces, each on a separate line
bcftools view -Ou -r NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11,NC_000024.10,NC_012920.1 $invcf | bcftools annotate -Oz -o $outvcf --rename-chrs $conversion -

