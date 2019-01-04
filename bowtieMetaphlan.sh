/home/praveen/Software/bowtie2-2.1.0/bowtie2 -p 8 -x mpa /home/praveen/Software/bowtie2-2.1.0/example/reads/ERR056989.fq > Experiment1.sam

/home/praveen/Software/bowtie2-2.1.0/bowtie2 --quiet --sam-no-hd --sam-no-sq --sensitive -x mpa -U /home/praveen/Software/bowtie2-2.1.0/example/reads/ERR056989.fq -S tmp.txt

./metaphlan.py --input_type blastout /home/praveen/Software/bowtie2-2.1.0/Experiment1.sam profiling_output.txt

./Software/metaphlan/metaphlan.py --bowtie2db bowtie2db/mpa --bt2_ps very-sensitive --input_type multifastq > profiled_samples/TD.txt

Software/metaphlan/metaphlan.py /home/praveen/Software/bowtie2-2.1.0/example/reads/ERR056989.fq -- bowtie2db mpa --nproc 5 --bowtie2out metagenome.bt2out.txt

Software/metaphlan/metaphlan.py home/praveen/Software/bowtie2-2.1.0/example/reads/ERR056989.fastq --bowtie2db bowtie2db/mpa --nproc 5 --bowtie2out metagenome.bt2out.txt

python Software/metaphlan/metaphlan.py --bowtie2db mpa home/praveen/Software/bowtie2-2.1.0/example/reads/ERR056989.fastq > 763577454-SRS014459-Stool.txt

Software/metaphlan/metaphlan.py --Software/metaphlan/bowtie2db mpa Software/bowtie2-2.1.0/example/reads/ERR056989.fastq > Stool.txt

Software/metaphlan/metaphlan.py --bowtie2db Software/metaphlan/bowtie2db/mpa Software/bowtie2-2.1.0/example/reads/ERR056989.fastq > Stool.txt

export BOWTIE2=/home/praveen/Software/bowtie2-2.1.0
export PATH=$BOWTIE2:$PATH
export PATH=$BOWTIE2/example/:$PATH
export PATH=$BOWTIE2/scripts/:$PATH
export PATH=$BOWTIE2/doc/:$PATH

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

mkdir output
utils/merge_metaphlan_tables.py profiled_samples/*.txt > output/merged_abundance_table.txt

mkdir output_images
plotting_scripts/metaphlan_hclust_heatmap.py -c bbcry --top 25 --minv 0.1 -s log --in output/merged_abundance_table.txt --out output_images/abundance_heatmap.png


mkdir tmp
metaphlan2graphlan.py profiled_samples/BF_1.txt  --tree_file tmp/BF_1.tree.txt --annot_file tmp/BF_1.annot.txt

~/Software/graphlan/graphlan_annotate.py --annot tmp/BF_1.annot.txt tmp/BF_1.tree.txt tmp/BF_1.xml

~/Software/graphlan/graphlan.py --dpi 200 tmp/BF_1.xml output_images/BF_1.png


$ mkdir -p tmp
$ for file in profiled_samples/*   
$ do
$     filename=`basename ${file}`
$     samplename=${filename%\.*}
$     plotting_scripts/metaphlan2graphlan.py ${file}  --tree_file tmp/${samplename}.tree.txt --annot_file tmp/${samplename}.annot.txt
$     graphlan_annotate.py --annot tmp/${samplename}.annot.txt tmp/${samplename}.tree.txt tmp/${samplename}.xml
$     graphlan.py --dpi 200 tmp/${samplename}.xml output_images/${samplename}.png
$ done



$ plotting_scripts/metaphlan2graphlan.py output/merged_abundance_table.txt  --tree_file tmp/merged.tree.txt --annot_file tmp/merged.annot.txt
$ graphlan_annotate.py --annot tmp/merged.annot.txt tmp/merged.tree.txt tmp/merged.xml
$ graphlan.py --dpi 200 tmp/merged.xml output_images/merged.png



$ sed 's/\([A-Z][A-Z]\)_\w*/\1/g' output/merged_abundance_table.txt > tmp/merged_abundance_table.4lefse.txt

The first LEfSe step consists of formatting the input table, making sure the class information is in the first row and scaling the values in [0,1M] which is useful for numerical computational reasons.

$ format_input.py tmp/merged_abundance_table.4lefse.txt tmp/merged_abundance_table.lefse -c 1 -o 1000000

Now, the LEfSe biomarker discovery tool can be used with default statistical options. Here we change one default parameter to increase the threshold on the LDA effect size from 2 (default) to 4.

$ run_lefse.py tmp/merged_abundance_table.lefse tmp/merged_abundance_table.lefse.out -l 4

The results of the operation can now be displayed plotting the resulting list of biomarkers with corresponding effect sizes.

$ plot_res.py --dpi 300 tmp/merged_abundance_table.lefse.out output_images/lefse_biomarkers.png

$ plot_cladogram.py --dpi 300 --format png tmp/merged_abundance_table.lefse.out output_images/lefse_biomarkers_cladogram.png 


$ plot_features.py -f one --feature_name "k__Bacteria.p__Firmicutes" tmp/merged_abundance_table.lefse tmp/merged_abundance_table.lefse.out output/Firmicutes.png

All features (or all biomarkers) can also be exported in one compressed archive:

$ plot_features.py -f diff --archive zip tmp/merged_abundance_table.lefse tmp/merged_abundance_table.lefse.out biomarkers.zip

$ mkdir output
$ utils/merge_metaphlan_tables.py humann_profiling/*.txt > output/merged_humann_abundance_table.txt




