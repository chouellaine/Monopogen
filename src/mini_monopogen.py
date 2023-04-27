import pysam
import os
import time
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


class MonopogenMini():
    """

    A class takes in a BAM file [bam_file] and a chromosome number of interest [chr] to run germline variant calling on the [chr] chunk of the BAM file. This is stripped down version of the Monopogen class

    Attributes
    ----------
    bam_file: [str] 
        file path of BAM file to run MonopogenMini
    out_fp: [str] 
        the file path to store output 
    chr_num: [int] 
        the chromosome number to run germline variant on
    imputation_panel_fp: [str] 
        the file path where the imputation panel is stored for the chromosome of interest, this is used to run Beagle
    reference_fp: [str] 
        the file path where the reference genome is stored
    app_path: [str]
        the file path where all dependent software executables are stored, these are bcftools, vcftools, bgzip, htslib, beagle.27Jul16.86a.jar. You may also use the folder "apps" of this repo
    
    Preconditions
    -------------
    the BAM file in bam_file must be sorted, prefixed with "Chr" annotation, and indexed
    
    Outputs
    -------
    * germline vcfs in the filepath defined in [out_fp]
    * log file

    """

    def __init__(self,
                bam_file,  
                out_fp,  
                chr_num,
                imputation_panel_fp,
                reference_fp,
                app_path):
        
        self.bam_file = bam_file
        self.chr_num = "chr{}".format(chr_num)
        self.out_fp = out_fp
        self.imputation_panel = imputation_panel_fp
        self.reference = reference_fp
        self.app_path = app_path
        self.max_mismatch = 4
        self.bam_filter = self.out_fp + "/Bam/" + self.chr_num + ".filter.bam"
        self.samtools = self.app_path + "/samtools"
        self.bcftools = self.app_path + "/bcftools"
        self.bgzip = self.app_path + "/bgzip"
        self.java = "java"
        self.beagle = self.app_path + "/beagle.27Jul16.86a.jar"
        return

    def BamFilter(self):
        logger.info("========== Calling BamFilter() ==========")
        start = time.time()
        infile = pysam.AlignmentFile(self.bam_file, "rb")
        tp = infile.header.to_dict()
        if not "RG" in tp:
            sampleID = os.path.splitext(os.path.basename(self.bam_file))[0]
            tp1 = [{'SM': sampleID, 'ID': sampleID,
                    'LB': "0.1", 'PL': "ILLUMINA", 'PU': sampleID}]
            tp.update({'RG': tp1})
        outfile = pysam.AlignmentFile(self.bam_filter, "wb", header=tp)
        for s in infile.fetch(self.chr_num):
            if s.has_tag("NM"):
                val = s.get_tag("NM")
            if s.has_tag("nM"):
                val = s.get_tag("nM")
            if val < self.max_mismatch:
                outfile.write(s)
        infile.close()
        outfile.close()
        os.system(self.samtools + " index " + self.out_fp + "/Bam/" + self.chr_num + ".filter.bam")
        end = time.time()
        time_elapsed = end - start
        logger.info("========== BamFilter Completion Time: {} \n".format(time_elapsed))
        self.bam_filter_time = time_elapsed
        return

    def SNVDiscover(self, filter_bam_fp=None):
        logger.info("========== Calling SNVDiscover() ==========")
        start = time.time()
        """
        HOTFIX: For some reason, the above command generates a vcf where the last chr position has the incorrect # of columns. Remove the last line of that vcf as a workaround.
        """
        cmd = self.samtools + " mpileup " + self.bam_filter + " -f " + \
            self.reference + " -r " + self.chr_num + " -q 20 -Q 20 -t DP -d 10000000 -v"
        cmd = cmd + " | " + self.bcftools + " view " + " | " + \
            self.bcftools + " norm -m-both -f " + self.reference
        cmd = cmd + " | grep -v \"<X>\" | grep -v INDEL | sed '$d' | " + \
            self.bgzip + " -c > " + self.out_fp + "/SCvarCall/" + self.chr_num + ".gl.vcf.gz"
        os.system(cmd)
        end = time.time()
        time_elapsed = end - start
        logger.info("========== SNVDiscover Completion Time: {} \n".format(time_elapsed))
        self.snv_discover_time = time_elapsed
        return

    def BeagleImpute(self, gl_fp=None):
        logger.info("========== Calling BeagleImpute() ==========")
        start = time.time()
        cmd = self.java + " -Xmx100g -jar " + self.beagle + " gl=" + self.out_fp + "/SCvarCall/" + self.chr_num + ".gl.vcf.gz" +  \
            " ref=" + self.imputation_panel + "  chrom=" + self.chr_num + " out=" + self.out_fp + "/SCvarCall/" + self.chr_num + ".gp " + \
            "impute=false  modelscale=2  nthreads=48  gprobs=true  niterations=0"
        os.system(cmd)
        end = time.time()
        time_elapsed = end - start
        logger.info("========== BeagleImpute Completion Time: {} \n".format(time_elapsed))
        self.impute_time = time_elapsed
        return

    def GenerateGermlineVCF(self):
        logger.info("========== Calling GenerateGermlineVCF() ==========")
        start = time.time()
        cmd = "zless -S " + self.out_fp + "/SCvarCall/" + self.chr_num + \
            ".gp.vcf.gz | grep -v  0/0  > " + self.out_fp + \
            "/SCvarCall/" + self.chr_num + ".germline.vcf"
        os.system(cmd)
        end = time.time()
        time_elapsed = end - start
        logger.info("========== GenerateGermlineVCF Completion Time: {} \n".format(
            time_elapsed))
        self.generate_germline_time = time_elapsed
        return

    def GenerateDirectory(self):
        out = self.out_fp
        os.system("mkdir -p " + out)
        os.system("mkdir -p " + out + "/Bam")
        os.system("mkdir -p " + out + "/SCvarCall")
        return

    def SCvarCall(self):
        logger.info("Running MonopogenMini\n")
        logger.info("BAM File: {}\n".format(self.bam_file))
        logger.info("CHR #: {}\n".format(self.chr_num))
        self.GenerateDirectory()
        self.BamFilter()
        self.SNVDiscover()
        self.BeagleImpute()
        self.GenerateGermlineVCF()
        return
