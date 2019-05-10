from pypedream.job import Job, required, optional, conditional


# FIXME: Could be adapted so it's possible to run directly from bam files, instead of using pre-calculated segments

class PureCN(Job):
    def __init__(self, input_seg=None, input_cnr=None, input_vcf=None, tumorid=None, outdir=None, funseg="none",
                 minpurity=0.05, hzdev=0.1, genome="hg19", postopt=False, maxnonclonal=None, minaf=0.02, error=0.001):
        Job.__init__(self)
        self.input_seg = input_seg
        self.input_cnr = input_cnr
        self.input_vcf = input_vcf
        self.tumorid = tumorid  # should be the same as the ID in the seg file
        self.outdir = outdir
        self.funseg = funseg
        self.minpurity = minpurity
        self.hzdev = hzdev
        self.genome = genome
        self.postopt = postopt
        self.maxnonclonal = maxnonclonal
        self.minaf = minaf
        self.error = error

        # FIXME: For PureCN jobs, the parameters *must* be specified to the constructor. This is inconsistent
        # across different Job sub-classes.

        # Determining output here purely as a means of determining when the job has
        # completed; PureCN itself determines the output file paths based on the
        # specified "--out" and "--sampleid" arguments.
        self.output = "{}/{}_loh.csv".format(
            self.outdir, self.tumorid)
        
        # Defining the names of the PureCN output files that need to be in place for downstream processing
        # These are not always produced by PureCN, and thus we need to "touch" the files after running PureCN
        self.out_csv = "{}/{}.csv".format(self.outdir, self.tumorid)
        self.out_genes = "{}/{}_genes.csv".format(self.outdir, self.tumorid)
        self.out_variants = "{}/{}_variants.csv".format(self.outdir, self.tumorid)
        # also _loh.csv needed, but it's defined above as output

        self.jobname = "purecn_{}".format(self.tumorid)

    def command(self):

        # activating conda env
        activate_cmd = "source activate purecn-env"

        # running PureCN
        running_cmd = "PureCN.R " + required("--out ", self.outdir) + \
                       required("--sampleid ", self.tumorid) + \
                       required("--segfile ", self.input_seg) + \
                       required("--tumor ", self.input_cnr) + \
                       required("--vcf ", self.input_vcf) + \
                       required("--genome ", self.genome) + \
                       optional("--funsegmentation ", self.funseg) + \
                       optional("--minpurity ", self.minpurity) + \
                       optional("--hzdev ", self.hzdev) + \
                       optional("--maxnonclonal ", self.maxnonclonal) + \
                       optional("--minaf ", self.minaf) + \
                       optional("--error ", self.error) + \
                       conditional(self.postopt, "--postoptimize")

        # deactivating the conda env
        deactivate_cmd = "conda deactivate"
        
        # touching required output files
        touch_cmd = "touch {} {} {} {}".format(self.out_csv, self.out_genes, self.out_variants, self.output)

        return " && ".join([activate_cmd, running_cmd, deactivate_cmd, touch_cmd])
