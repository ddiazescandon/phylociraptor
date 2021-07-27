configfile: "data/config.yaml"

import os
import glob

def get_assemblies(wildcards):
	sp = "{wildcards.species}".format(wildcards=wildcards)
	sp.replace(" ", "_")
	#print(sp)
	if os.path.isfile("results/assemblies/" + sp + ".fna"):
		return ["results/assemblies/" + sp + ".fna"]
	elif os.path.isfile("results/assemblies/" + sp + ".fna.gz"):
		return ["results/assemblies/" + sp + ".fna.gz"]
	else:
		return []

rule run_busco:
	input:
		#assembly = "results/assemblies/{species}.fna",
		assembly = get_assemblies,
		augustus_config_path = "results/augustus_config_path",
		busco_set = "results/busco_set/"+config["busco"]["set"]
	output:
		checkpoint = "results/checkpoints/busco/busco_{species}.done",
		augustus_output = "results/busco/{species}/run_busco/augustus_output.tar.gz",
		blast_output = "results/busco/{species}/run_busco/blast_output.tar.gz",
		hmmer_output = "results/busco/{species}/run_busco/hmmer_output.tar.gz",
		logs = "results/busco/{species}/run_busco/logs.tar.gz",
		# augustus_model = "results/busco/{species}/run_busco/{species}_augustus_model.tar.gz",
		full_table = "results/busco/{species}/run_busco/full_table_busco.tsv",
		short_summary ="results/busco/{species}/run_busco/short_summary_busco.txt",
		missing_busco_list ="results/busco/{species}/run_busco/missing_busco_list_busco.tsv",
		single_copy_buscos = "results/busco/{species}/run_busco/single_copy_busco_sequences.tar",
		single_copy_buscos_tarlist = "results/busco/{species}/run_busco/single_copy_busco_sequences.txt"

	benchmark: "results/statistics/benchmarks/busco/run_busco_{species}.txt"
	threads: int(config["busco"]["threads"])
	shadow: "shallow"
	log: "log/{species}_busco.log"
	params:
		wd = os.getcwd(),
		sp = config["busco"]["augustus_species"],
		additional_params = config["busco"]["additional_parameters"],
		species = lambda wildcards: wildcards.species,
		mode = "genome",
		augustus_config_in_container = "/usr/local/config",
		set = config["busco"]["set"]
	singularity:
		"docker://ezlabgva/busco:v5.2.1_cv1"
	shell:
		"""
		mkdir -p log
		dir=results/busco/{params.species}
		# prepare stripped down version auf augustus config path.
		# this is introduced to lower the number of files.
		mkdir augustus
		cp -R {params.augustus_config_in_container}/cgp augustus
		cp -R {params.augustus_config_in_container}/extrinsic augustus
		cp -R {params.augustus_config_in_container}/model augustus
		cp -R {params.augustus_config_in_container}/profile augustus
		mkdir augustus/species
		cp -R {params.augustus_config_in_container}/species/generic augustus/species/
		
		if [ -d {params.augustus_config_in_container}/species/{params.sp} ]
		then
			cp -R {params.augustus_config_in_container}/species/{params.sp} augustus/species
		fi		

		export AUGUSTUS_CONFIG_PATH=$(pwd)/augustus
		#export AUGUSTUS_SCRIPTS_PATH=/usr/local/bin/ #might be necessary in ezlabgva/busco:v5.2.1_cv1 image
		#export AUGUSTUS_BIN_PATH=/usr/local/bin/
		echo $AUGUSTUS_CONFIG_PATH
		
		# handle gzipped and other assemblies differently
		if [[ "{input.assembly}" =~ \.gz$ ]]
		then
			fullname=$(basename "{input.assembly}")
			filename="${{fullname%.*}}"
			gunzip -c $(readlink -f "{input.assembly}") > "$filename"
		else
			filename="{input.assembly}"
		fi
		#echo "Assembly used for BUSCO is "$filename 2>&1 | tee {log}
		busco -i $filename -f --out {params.species} -c {threads} --augustus --augustus_species {params.sp} --lineage_dataset $(pwd)/{input.busco_set} -m {params.mode} {params.additional_params} 2>&1 | tee {log}
		# do some cleanup to save space
		basedir=$(pwd)
		cd {params.species}/run_{params.set}
		$basedir/bin/tar_folder.sh $basedir/{output.blast_output} blast_output
		$basedir/bin/tar_folder.sh $basedir/{output.hmmer_output} hmmer_output
		$basedir/bin/tar_folder.sh $basedir/{output.augustus_output} augustus_output
		cd ..
		$basedir/bin/tar_folder.sh $basedir/{output.logs} logs
		cd ..
		tar -pcf {output.single_copy_buscos} -C {params.species}/run_{params.set}/busco_sequences single_copy_busco_sequences 
		tar -tvf {output.single_copy_buscos} > {output.single_copy_buscos_tarlist}	

		#move output files:
		mv {params.species}/run_{params.set}/full_table.tsv {output.full_table}
		mv {params.species}/run_{params.set}/short_summary.txt {output.short_summary}
		mv {params.species}/run_{params.set}/missing_busco_list.tsv {output.missing_busco_list}
		
		buscos=$(tail -n +6 {output.full_table} | cut -f 2 | sort | uniq -c | tr '\\n' ' ' | sed 's/ $/\\n/g')
		name="{params.species}"
		echo "$(date) $name $buscos" >> results/statistics/runlog.txt
		
		#touch checkpoint
		touch {output.checkpoint}
		"""

rule busco:
	input:
		checks = expand("results/checkpoints/busco/busco_{species}.done", species=glob_wildcards("results/assemblies/{species}.fna").species + glob_wildcards("results/assemblies/{species}.fna.gz").species)
	output:
		"results/checkpoints/busco.done"
	shell:
		"""
		touch {output}
		"""

rule extract_busco_table:
	input:
		busco_set = "results/busco_set/"+config["busco"]["set"],
		busco = rules.busco.output
	output:
		busco_table = "results/busco_table/busco_table.txt",
		#checkpoint = "results/checkpoints/extract_busco_table.done"
	benchmark:
		"results/statistics/benchmarks/busco/extract_busco_table.txt"
	singularity:
		"docker://reslp/biopython_plus:1.77"
	params:
		busco_dir = "results/busco/"
	shell:
		"""
		python bin/extract_busco_table.py --hmm {input.busco_set}/hmms --busco_results {params.busco_dir} -o {output.busco_table}
		echo "species\tcomplete\tsingle_copy\tduplicated\tfragmented\tmissing\ttotal" > results/statistics/busco_summary.txt
		for file in $(ls results/busco/*/run_busco/short_summary_busco.txt);do  name=$(echo $file | sed 's#results/busco/##' | sed 's#/run_busco/short_summary_busco.txt##'); printf $name; cat $file | grep -P '\t\d' | awk -F "\t" '{{printf "\t"$2}}' | awk '{{print}}'; done >> results/statistics/busco_summary.txt
		"""
rule orthology:
	input:
		expand("results/checkpoints/busco/busco_{species}.done", species=glob_wildcards("results/assemblies/{species}.fna").species),
		#"results/checkpoints/busco.done",
		"results/busco_table/busco_table.txt",
		#"results/checkpoints/create_sequence_files.done",
		#"results/checkpoints/remove_duplicated_sequence_files.done"
	output:
		"checkpoints/orthology.done"
	shell:
		"""
		touch {output}
		echo "$(date) - Pipeline part 1 (orthology) done." >> results/statistics/runlog.txt
		"""
