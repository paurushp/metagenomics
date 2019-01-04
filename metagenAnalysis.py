#####################################################################################
## Description: Python code to analyse Metagenomic data				    #
## Author: Paurush Praveen							    #
## Contact: praveen@cosbi.eu							    #
## Biology: 16s sequences from infant infant stool samples with two 	  	    #
## diets 1. Breast Fed and 2. Formula Fed. Formula fed assumed as control	    #
#####################################################################################

# Install and add metaphlan+BOWTIE2 to PATH

export BOWTIE2=/home/praveen/Software/bowtie2-2.1.0
export PATH=$BOWTIE2:$PATH
export PATH=$BOWTIE2/example/:$PATH
export PATH=$BOWTIE2/scripts/:$PATH
export PATH=$BOWTIE2/doc/:$PATH

# ## Run metaphlan on all the file 

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056989.fastq > Stool1.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056990.fastq > Stool2.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056991.fastq > Stool3.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056992.fastq > Stool4.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056993.fastq > Stool5.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056994.fastq > Stool6.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056995.fastq > Stool7.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056996.fastq > Stool8.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056997.fastq > Stool9.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056998.fastq > Stool10.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056999.fastq > Stool11.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR057000.fastq > Stool12.txt


# ### Merge the results files ***** USE IT FOR MERGING GENE NAME-MICROBE  RELATIONS AS WELL *****
mkdir output
utils/merge_metaphlan_tables.py profiled_samples/*.txt > output/merged_abundance_table.txt

# ### Plot the heatmap for the merged file

mkdir output_images
plotting_scripts/metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log --in output/merged_abundance_table.txt --out output_images/abundance_heatmap.png


# ### All files in a directory
mkdir -p tmp
 for file in profiled_samples/*   
 do
     filename=`basename ${file}`
     samplename=${filename%\.*}
     plotting_scripts/metaphlan2graphlan.py ${file}  --tree_file tmp/${samplename}.tree.txt --annot_file tmp/${samplename}.annot.txt
     graphlan_annotate.py --annot tmp/${samplename}.annot.txt tmp/${samplename}.tree.txt tmp/${samplename}.xml
     graphlan.py --dpi 200 tmp/${samplename}.xml output_images/${samplename}.png
done


# ### For the merged annotations
 plotting_scripts/metaphlan2graphlan.py output/merged_abundance_table.txt  --tree_file tmp/merged.tree.txt --annot_file tmp/merged.annot.txt
 ~/Software/graphlan/graphlan_annotate.py --annot tmp/merged.annot.txt tmp/merged.tree.txt tmp/merged.xml
 ~/Software/graphlan/graphlan.py --dpi 200 tmp/merged.xml output_images/merged.png

##############################################################
plotting_scripts/metaphlan2graphlan.py output/merged_abundance_tableBF.txt  --tree_file tmp/merged.treeBF.txt --annot_file tmp/merged.annotBF.txt
~/Software/metaphlan$  ~/Software/graphlan/graphlan_annotate.py --annot tmp/merged.annotBF.txt tmp/merged.treeBF.txt tmp/mergedBF.xml
~/Software/metaphlan$  ~/Software/graphlan/graphlan.py --dpi 200 tmp/mergedBF.xml output_images/mergedBF.png


plotting_scripts/metaphlan2graphlan.py output/merged_abundance_tableFF.txt  --tree_file tmp/merged.treeFF.txt --annot_file tmp/merged.annotFF.txt
 ~/Software/graphlan/graphlan_annotate.py --annot tmp/merged.annotFF.txt tmp/merged.treeFF.txt tmp/mergedFF.xml
 ~/Software/graphlan/graphlan.py --dpi 200 tmp/mergedFF.xml output_images/mergedFF.png


