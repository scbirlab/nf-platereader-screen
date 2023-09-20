#!/usr/bin/env nextflow

/*
========================================================================================
   Platereader screen Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-platereader-screen
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println """\
         S C B I R   P L A T E R E A D E R   S C R E E N   P I P E L I N E
         =================================================================

         Nextflow pipeline to process files from Biotek platereaders, normalise 
         output, perform QC, make plots, and do statistical tests. 

         Usage:
            nextflow run sbcirlab/nf-platereader-screen --sample_sheet <csv> --data_dir <dir> --output_dir <dir> --export_layout (plate|row) --positive <neg> --negative <pos> --grouping <cols> --hit_grouping <cols>
            nextflow run sbcirlab/nf-platereader-screen -c <config-file>

         Required parameters:
            sample_sheet      Path to a CSV containing sample IDs matched with ONT barcode IDs, genome information, and adapter sequences.
            data_dir          Path to directory where files in `sample_sheet` can be found.
            output_dir        Path to directory where outputs should be saved.
            export_layout     "plate" or "row", indicating the layout of the exported platereader data.
            positive          Positive control label contained in the `control_column`.
            negative          Negative control label contained in the `control_column`.
            grouping          Columns with labels indicating groups (i.e batches) within which to normalize data. Often this will be a column of plate labels.
            hit_grouping      Columns with labels indicating repeated measurements of experimental conditions. Often this will be a column indicating different compounds.

         Optional parameters (with defaults):   
            layout_content_name = "compound_name"    A name to give the data from the wells of the compound source plates.
            control_column = "compound_name"         Name of column containing labels indicating positive and negative controls.

         Optional parameters:
            counterscreen     A colon-separated column name and negative control label to do a second wave of hit calling. For example `"strain_name:WT"`.

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}
if ( !params.data_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to data_dir")
}
if ( !params.output_dir ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to output_dir")
}
if ( !params.export_layout ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide the export_layout: 'row' or 'plate'")
}
if ( !params.positive ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a positive control label.")
}
if ( !params.negative ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a negative control label.")
}
if ( !params.grouping ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a batch grouping.")
}
if ( !params.hit_grouping ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a hit_grouping (groups of experimental conditions to do statstical tests).")
}


wd = file(params.sample_sheet)

normalized_o = "${params.output_dir}/normalized"
qc_o = "${params.output_dir}/qc"
plots_o = "${params.output_dir}/plots"
hits_o = "${params.output_dir}/hits"

log.info """\
         S C B I R   P L A T E R E A D E R   S C R E E N   P I P E L I N E
         =================================================================
         inputs
            sample sheet   : ${params.sample_sheet}
            input dir      : ${params.data_dir}
         data
            export layout  : ${params.export_layout}
         layout  
            content name   : ${params.layout_content_name}
         controls
            control column : ${params.control_column}
            positive       : ${params.positive}
            negative       : ${params.negative}
            grouping       : ${params.grouping}
         hit calling 
            grouping       : ${params.hit_grouping}
            counterscreen  : ${params.counterscreen}
         output            : ${params.output_dir}
            normalized     : ${normalized_o}
            qc             : ${qc_o}
            hits           : ${hits_o}
            plots          : ${plots_o}
         """
         .stripIndent()


dirs_to_make = [normalized_o, qc_o, hits_o, plots_o]
log.info  """
          Making directories:
          """.stripIndent()

dirs_to_make.each { 
   log.info "${it}: " 
   log.info file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}
/*
========================================================================================
   Create Channels
========================================================================================
*/

sample_sheet_ch = Channel.fromPath( params.sample_sheet, 
                                    checkIfExists: true )

sample_ch = sample_sheet_ch
                   .splitCsv( header: true )
                   .map { tuple(it.experiment_id, 
                                file("${params.data_dir}/${it.data_filename}")) }
                   .groupTuple()
                   .map { key, words -> tuple( groupKey(key, words.size()), words ) }

layout_ch = sample_sheet_ch
                   .splitCsv( header: true )
                   .map { tuple(it.experiment_id,
                                file("${params.data_dir}/${it.compound_source_filename}")) }
                   .unique()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   DATA2COLUMNS_and_NORMALIZE(sample_ch.join(layout_ch).combine(sample_sheet_ch))

   DATA2COLUMNS_and_NORMALIZE.out | QC
   DATA2COLUMNS_and_NORMALIZE.out | PLOTS
   DATA2COLUMNS_and_NORMALIZE.out | CALL_HITS

   if ( params.counterscreen ) {
      DATA2COLUMNS_and_NORMALIZE.out | CALL_HITS_COUNTERSCREEN 
   } 

}

/*
 * Convert the Gen5 exported file to columns.
 */
process DATA2COLUMNS_and_NORMALIZE {

   tag "${expt_id}" 

   publishDir( normalized_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data ), path( compound_source_layout ), path( sample_sheet )

   output:
   tuple val( expt_id ), path( "*_normalized.tsv" )

   script:
   """
   hts pivot ${compound_source_layout} \
      --name ${params.layout_content_name} \
      --prefix compound_source \
      > ${expt_id}-compounds.tsv

   hts parse ${data} --data-shape ${params.export_layout} \
      | hts join --right ${sample_sheet} \
      | hts join --right ${expt_id}-compounds.tsv \
      | hts normalize \
         --control ${params.control_column} \
         --positive ${params.positive} \
         --negative ${params.negative} \
         --grouping ${params.grouping} \
         --output ${expt_id}_normalized.tsv
   """

}

/*
 * Merge sample sheet with data files.
 */
process QC {

   tag "${expt_id}"

   publishDir( qc_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data )

   output:
   tuple val( expt_id ), path( "*.tsv" ), path( "*.png" )

   script:
   """
   hts qc ${data} \
      --control ${params.control_column} \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --grouping ${params.grouping} \
      --plot ${expt_id} \
      --output ${expt_id}_qc.tsv
   """

}

/*
 * Summarize replicates with statistics.
 */
process CALL_HITS {

   tag "${expt_id}"

   publishDir( hits_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data )

   output:
   tuple val( expt_id ), path( "*.tsv" ), path( "*.png" )

   script:
   """
   hts summarize ${data} \
      --control ${params.control_column} \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --grouping ${params.hit_grouping} \
      --plot ${expt_id}_summary \
      --output "${expt_id}_summary.tsv" 
   """

}

/*
 * Summarize replicates with statistics between data subgroups.
 */
process CALL_HITS_COUNTERSCREEN {

   tag "${expt_id}"

   publishDir( hits_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data )

   output:
   tuple val( expt_id ), path( "*.tsv" ), path( "*.png" )

   script:
   """
   hts summarize ${data} \
      --control ${params.counterscreen.split(":")[0]} \
      --positive ${params.positive} \
      --negative ${params.counterscreen.split(":")[1]} \
      --grouping ${params.hit_grouping} \
      --plot ${expt_id}_summary-counter \
      --output "${expt_id}_summary-counter.tsv" 
   """

}

/*
 * Visualize plate data.
 */
process PLOTS {

   tag "${expt_id}"

   publishDir( plots_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data )

   output:
   tuple val( expt_id ), path( "*.png" )

   script:
   """
   hts plot-hm ${data} \
      --grouping ${params.grouping} \
      --output "${expt_id}" 

   hts plot-rep ${data} \
      --control ${params.control_column} \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --grouping ${params.hit_grouping} \
      --output "${expt_id}" 

   hts plot-hist ${data} \
      --control ${params.control_column} \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --output "${expt_id}" 

   """

}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/