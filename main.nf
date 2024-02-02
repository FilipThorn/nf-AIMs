#!/usr/bin/env nextflow

if (params.help) {
    log.info """\
            __________________
            |                |
            | |```````| |`````
            | |____   | |
            |     |   | |
            | |````   | |
            | |ilip   | |hörn     
            –––––––––––––––––––––––––––––––––––––––
            AIMs  
            NEXTFLOW   P I P E L I N E                
            –––––––––––––––––––––––––––––––––––––––
            'USAGE'
            nextflow run main.nf 
         
            'Mandatory arguments:'
            --outdir       PATH      Path to output directory where results will be should be stored 
            
            'OPTIONS'
            --help                   Outputs this help log      
            -resume                  Nextflow cmd to resume modified workflow
            
            'HPC'
            -profile       FILE      If intention to run workflow on HPC please provide a suitable profile 
                                     in the nextflow.config file 

            'SUPPORT'
            Email Filip.Thorn@NRM.se for questions on script
            """
    exit 1
}


log.info """\
         –––––––––––––––––––––––––––––––––––––––
         Ancestry Painter 
         NEXTFLOW   P I P E L I N E                
         –––––––––––––––––––––––––––––––––––––––
         outdir       	: ${params.outdir}
         """
         .stripIndent()


// Channel
	channel.fromPath(params.hybrid_tsv)
		.splitCsv(header:true, sep:'\t')
		.view()
		.map { row -> tuple(row.Hybrid, file(row.VCF), file(row.Parent1), file(row.Parent2) ) }
		.set { input_ch }


process CalcFst {
   
   tag "$hybrid"   
    
   publishDir "${params.outdir}/01.Fst/$hybrid", mode:'copy'

   input:
   tuple val(hybrid), file(vcf), file(parent1), file(parent2) from input_ch

   output:
   tuple val(hybrid), file("${hybrid}_${parent1.baseName}.weir.fst"), file("${hybrid}_${parent2.baseName}.weir.fst"), file("${parent1.baseName}_${parent2.baseName}.weir.fst")  into fst_ch

   script:
   """
    echo $hybrid > hyb.txt

    vcftools --gzvcf ${vcf} --weir-fst-pop ${parent1}  --weir-fst-pop ${parent2} --out ${parent1.baseName}_${parent2.baseName}
    vcftools --gzvcf ${vcf} --weir-fst-pop hyb.txt  --weir-fst-pop ${parent1} --out ${hybrid}_${parent1.baseName}
    vcftools --gzvcf ${vcf} --weir-fst-pop hyb.txt  --weir-fst-pop ${parent2} --out ${hybrid}_${parent2.baseName}

   """
}


process CalcAIMs {

   label 'RAM_high'    

   tag "$hybrid"
   
   publishDir "${params.outdir}/02.AIMs/$hybrid", mode:'copy'

   input:
   tuple val(hybrid), file(HP1), file(HP2), file(P1P2) from fst_ch

   output:
   tuple val(hybrid), file("${hybrid}_AIMs.csv"), val("${P1P2.baseName.tokenize(".")[0].tokenize("_")[0]}"), val("${P1P2.baseName.tokenize(".")[0].tokenize("_")[1]}") into AIMs_ch

   script:
   """
   python $params.AIMs_path --Parents_fst ${P1P2} --P1H_fst ${HP1} --P2H_fst ${HP2} --P1 ${HP1.baseName.tokenize(".")[0].tokenize("_")[-1]} --P2 ${HP2.baseName.tokenize(".")[0].tokenize("_")[-1]} --Out ${hybrid}_AIMs.csv

   """
}


process PlotAIMs {

   label 'FAST'

   tag "$hybrid"   

   publishDir "${params.outdir}/03.plot/", mode:'copy'

   input:
   tuple val(hybrid), file(csv), val(P1), val(P2) from AIMs_ch

   output:
   file("*")

   script:
   """




   #plot with r
   Rscript $params.plot_path $csv $hybrid $P1 $P2  

   """
}

