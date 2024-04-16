#!/usr/bin/env nextflow
/*
Copyright 2024 Francisco Falcon (francisco .falcon [at] imp.ac.at)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and
or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
*/
nextflow.enable.dsl=2
params.Reads      =""
params.CellTagSeq ="GGGTACCGCGGGCCCG"
params.TSOReadSeq ="TTTCTTATATGGG"
params.chunkSize  =100000
params.V1Barcode  ="ACCGGT"
params.V2Barcode  ="CATCAC"
params.V3Barcode  ="CGTACA"
params.whitelist  =""
params.help       =""
params.outDir     =""
def helpMessage() {
  log.info"""
  =======================================================================
    ╔═╗┌─┐┬  ┬ ╔╦╗┌─┐┌─┐    ╔╗ ┌─┐┬─┐┌─┐┌─┐┌┬┐┌─┐┌─┐
    ║  ├┤ │  │  ║ ├─┤│ ┬    ╠╩╗├─┤├┬┘│  │ │ ││├┤ └─┐
    ╚═╝└─┘┴─┘┴─┘╩ ┴ ┴└─┘────╚═╝┴ ┴┴└─└─┘└─┘─┴┘└─┘└─┘
                                                    by Francisco Falcon®
  =======================================================================
  Usage: nextflow run nf-TE_RNASeq.nf [PARAMETERS]
  Parameters:
  --Reads          Read file in fq.gz format
  --TSOReadSeq     Sequence upstream of the 10x barcode (Default:GGGTACCGCGGGCCCG )
  --CellTagSeq     Sequence upstream of the CellTag Barcode (Default:TTTCTTATATGGG )
  --chunkSize      Read size to split the fastq file into (Default:5000)
  --V1Barcode      Barcode sequence of the V1 barcode (Default:ACCGGT)
  --V2Barcode      Barcode sequence of the V2 barcode (Default:CATCAC)
  --V3Barcode      Barcode sequence of the V3 barcode (Default:CGTACA)
  --outDir         Output directory (Default:Current Directory)
  --whitelist      10x whitelist file
  --help           This beautiful help message :)
  """.stripIndent()
  return ""
}
if (params.help){
    exit 0, helpMessage()
}
if ( !params.Reads){
    exit 1, helpMessage() + "--Reads not specified !!!"
} else if (!params.whitelist){
    exit 1, helpMessage() + "--whitelist not specified !!!"
}
log.info " Running the pipeline with the following parameters: "
log.info "====================================================="
log.info "Read File      : ${params.Reads}"
log.info "TSOReadSeq     : ${params.TSOReadSeq}"
log.info "CellTagSeq     : ${params.CellTagSeq}"
log.info "chunkSize      : ${params.chunkSize}"
log.info "V1Barcode      : ${params.V1Barcode}"
log.info "V2Barcode      : ${params.V2Barcode}"
log.info "V3Barcode      : ${params.V3Barcode}"
log.info "10x whitelist  : ${params.whitelist}"
//Get absolute paths
Reads    =file(params.Reads)
Whitelist=file(params.whitelist)
outDir   =file(params.outDir)

//Split reads into chunks
process ReadSplit{
    executor = 'slurm'
    cpus = 8
    time = { 2.day }
    memory = { 10.GB }
    clusterOptions = {'--qos=medium --partition=c'}
    output:
    path("*.fa.gz"), emit: ReadChunksFA
    path("*.fastq.gz"), emit: ReadChunksFQ
    script:
    """
    #Split Files
    seqkit split2 --threads ${task.cpus} -s ${params.chunkSize} -e ".gz" -1 ${Reads} -O ./ --by-size-prefix "SplitReads"
    #Convert to Fasta
    for FqFiles in \$(ls *fastq.gz);do
      BaseName=\${FqFiles%%.*}
      seqkit fq2fa \${FqFiles} --threads ${task.cpus} | seqkit replace -p ":" -r "__"| seqkit replace -p "\s.+" > \${BaseName}.fa.gz
    done
    """
}

//Search for 10x Barcodes
process Barcode_10xSearch{
    executor = 'slurm'
    cpus = 1
    time = { 1.hour }
    memory = { 10.GB }
    errorStrategy = { 'retry' }
    maxRetries = 3
    clusterOptions = {'--qos=rapid --partition=c'}
    input:
    tuple val(ReadBaseName),path(ReadsFQ),path(ReadsFA)
    output:
    tuple val(ReadBaseName),path("BarcodeLocation.TenX.txt"),path("BarcodeSeq.TenX.txt")
    script:
    """
    #Locate Barcode in Reads
    zcat ${ReadsFQ} | seqkit locate -m 1 -i -p ${params.TSOReadSeq} > BarcodeLocation.TenX.txt
    #Locate Barcode
    samtools faidx -r <(cat BarcodeLocation.TenX.txt | sed "1d" | cut -f1,5,6| awk '{OFS="\\t"}\$2>= 27{\$3=\$2-1;\$2=\$3-25;print}' |sed 's/:/__/g;s/\\t/:/;s/\\t/-/') ${ReadsFA}| seqkit fx2tab| awk '{OFS="\\t"}{print \$1,substr(\$2,1,16),substr(\$2,17,10)}' > BarcodeSeq.TenX.txt
    """
}

//Search for CellTag Barcodes
process Barcode_CellTagSearch{
    executor = 'slurm'
    cpus = 1
    time = { 1.hour }
    memory = { 10.GB }
    errorStrategy = { 'retry' }
    maxRetries = 3
    clusterOptions = {'--qos=rapid --partition=c'}
    input:
    tuple val(ReadBaseName),path(ReadsFQ),path(ReadsFA)
    output:
    tuple val(ReadBaseName),path("BarcodeLocation.CellTag.txt"),path("BarcodeSeq.CellTag.txt")
    """
    #Locate Barcode in Read
    zcat ${ReadsFQ} | seqkit locate -m 1 -i -P -p ${params.CellTagSeq} | awk '\$6 < 135{print}' > BarcodeLocation.CellTag.txt
    
    #Extract barcode sequence
    samtools faidx -r <(cat BarcodeLocation.CellTag.txt | cut -f1,5,6| awk '{OFS="\\t"}{\$2=\$3+1;\$3=\$2+14;print}'| sed 's/:/__/g;s/\\t/:/;s/\\t/-/') ${ReadsFA} | seqkit fx2tab| awk '{OFS="\\t"}{print \$1,substr(\$2,1,9),substr(\$2,10,6)}' > BarcodeSeq.CellTag.txt
    """
}
//Combine barcodes
process Barcode_Joining{
    executor = 'slurm'
    cpus = 1
    time = { 1.hour }
    memory = { 5.GB }
    errorStrategy = { 'retry' }
    maxRetries = 3
    clusterOptions = {'--qos=rapid --partition=c'}

    input:
    tuple val(ReadBaseName),path(BarcodeTenX_Loc),path(BarcodeTenX_Seq),path(BarcodeCellTag_Loc),path(BarcodeCellTag_Seq)
    output:
    path('BarcodeTable.CellTag.10X.tsv')
    script:
    """
    gawk '{OFS="\\t"}
    {Length=\$1
    gsub(".*:","",Length)
    gsub("-.*","",Length)
    gsub(":.*","",\$1)
    if(NR==FNR){
        Matrix[\$1][1]=\$2
        Matrix[\$1][2]=\$3
        Matrix[\$1][3]=Length
        Matrix[\$1][4]=Length+length(\$2)-1
        Matrix[\$1][5]=Length+length(\$2)
        Matrix[\$1][6]=Length+length(\$2)+length(\$3)-1
        for (i=7; i<=12;i++) Matrix[\$1][i]="."
    }else{
        Matrix[\$1][7]=\$2
        Matrix[\$1][8]=\$3
        Matrix[\$1][9]=Length
        Matrix[\$1][10]=Length+length(\$2)-1
        Matrix[\$1][11]=Length+length(\$2)
        Matrix[\$1][12]=Length+length(\$2)+length(\$3)-1
        if(Matrix[\$1][1]==""){
            for (i=1; i<=6;i++) Matrix[\$1][i]="."
        }
    }
    }END{for (key in Matrix){
        printf key"\\t"
        for (i=1;i<=12;i++) printf Matrix[key][i]"\\t"
        print ""
    }
    }' ${BarcodeCellTag_Seq} ${BarcodeTenX_Seq} > BarcodeTable.tsv
    #Select only for Barcode Versions
    gawk '{OFS="\\t"}
    \$3=="${params.V1Barcode}" || \$3 =="${params.V2Barcode}" || \$3=="${params.V3Barcode}"{
        print \$0
        }' BarcodeTable.tsv > BarcodeTable.CellTag.tsv
    #Select 10x Barcodes that are within the whitelist
    gawk '{OFS="\\t"}
    {if(NR==FNR){
        Whitelist[\$1]=\$1
    } else{
        if (\$8 in Whitelist){
            print \$0
        }
    }
    }' ${Whitelist} BarcodeTable.CellTag.tsv | sed '1i ReadID\\tCellTagBarcode\\tCellTagVersion\\tCTBarcode_Start\\tCTBarcode_End\\tCTVersion_Start\\tCTVersion_End\\tCellID\\tUMI\\tCellIDStart\\tCellIDEnd\\tUMIStart\\tUMIEnd'> BarcodeTable.CellTag.10X.tsv

    #Extract only the columns containing barcodes and UMIs, remove repeated entries
    cut -f2,3,8,9 BarcodeTable.CellTag.10X.tsv | sort | uniq -c > Temp.txt

    #Filter by counts 
    cut -f1-3 temp.txt| sort | uniq -c | awk '{OFS="\t"}$1>=2{print $1,$2,$3,$4}' > Barcode.Uniq.FilteredCounts.txt
    """
}

//Workflow
workflow {
    ReadSplit()
    FastaFiles=ReadSplit.out.ReadChunksFA.flatten().map{n -> [n.simpleName,n]}
    FastqFiles=ReadSplit.out.ReadChunksFQ.flatten().map{n -> [n.simpleName,n]}
    ReadFiles=FastqFiles.join(FastaFiles)
    Barcode_CellTagSearch(ReadFiles)
    Barcode_10xSearch(ReadFiles)
    CombinedBarcodes=Barcode_10xSearch.out.join(Barcode_CellTagSearch.out)
    Barcode_Joining(CombinedBarcodes)
    Barcode_Joining.out.collectFile(keepHeader:true,storeDir:"${outDir}",name:"BarcodesOut.txt")
}
