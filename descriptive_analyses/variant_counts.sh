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
