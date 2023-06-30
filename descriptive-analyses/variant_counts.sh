# 2023-06-16
# Gwenna Breton
# Goal: count known and unreported variants in our newly generated data using an updated version of dbsnp.

# Preliminary step: retrieve and prepare dbsnp file
ssh gwennabr@rackham.uppmax.uu.se
screen -S dbsnp
cd /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156
interactive -A p2018003 -p core -n 2 -t 2:00:00 -M snowy
##To list the files:
curl ftp://ftp.ncbi.nih.gov/snp/redesign/archive/b156/VCF/
##To download the files:
wget ftp://ftp.ncbi.nih.gov/snp/redesign/archive/b156/VCF/GCF_000001405.40.gz.md5 #And all the other files in the folder

#The relevant file for hg38 is GCF_000001405.40.gz. However the contigs are named like NC_000001.11 for chr1, etc. I need to change the contig names to what I have in my VCF files.
# I decided to keep only 1-22XYMT (I could skip Y and MT actually) and rename these to what I have in my files.
cd /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/modified_version
zgrep -v "##" ../VCF/GCF_000001405.40.gz | cut -f1 | sort | uniq > list_of_contigs
grep NC list_of_contigs > list_of_contigs_conversion_1-22XYMT
#Then I made a list that listed the "new" names and pasted the two lists together.

##Select 1-22XYMT and rename the contigs (trying to do all in one go)

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

# Count variants in the VCF files after VQSR, relatedness, HWE and geno0.1 filtering
#This is based on code in /Users/gwennabreton/Documents/Previous_work/PhD_work/P2_RHG_KS/sequence_processing/scripts/mapping_and_GATK_final_scripts/HC_BPresolution/initial_analyses.sh
cd /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/bash_outputs/2023
#chrom=22 #Test run with chr22
for chrom in {1..21}; do
(echo '#!/bin/bash -l'
echo "
module load bioinfo-tools picard/2.10.3
prefix=25KS.48RHG.104comp.HCBP.${chrom}.recalSNP99.9.recalINDEL99.0.FAIL1FAIL2FAIL3.reheaded
folder=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/HC_BPresolution/3maskrecal.realn/allsites/3_geno01_hwefiltering/
cd \$SNIC_TMP
cp \${folder}\${prefix}.vcf.gz .
cp \${folder}\${prefix}.vcf.gz.tbi .
cp /crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/dbsnp156/modified_version/dbsnp156_1-22XYMT.vcf.gz .
java -Xmx6g -jar \$PICARD_HOME/picard.jar CollectVariantCallingMetrics \
INPUT=\${prefix}.vcf.gz \
OUTPUT=\${prefix}.dbsnp156 \
DBSNP=dbsnp156_1-22XYMT.vcf.gz
cp \${prefix}.dbsnp156* \${folder}
exit 0") | sbatch -p core -n 1 -t 36:0:0 -A p2018003  -J filtered_${chrom}_CVCM -o filtered_${chrom}_CVCM_dbsnp156.output -e filtered_${chrom}_CVCM_dbsnp156.output --mail-user gwenna.breton@ebc.uu.se --mail-type=END,FAIL -M snowy
done
#To see stats for snowy: jobinfo -u gwennabr -M snowy
