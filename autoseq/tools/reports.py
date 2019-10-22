import os
import uuid

from pyjobmanager.job import Job, required, conditional


class CompileMetadata(Job):
    def __init__(self, referral_db_conf, blood_barcode, tumor_barcode, output_json, addresses):
        Job.__init__(self)
        self.referral_db_conf = referral_db_conf
        self.blood_barcode = blood_barcode
        self.tumor_barcode = tumor_barcode
        self.output_json = output_json
        self.addresses = addresses

    def command(self):
        # compileMetadata 3098121 3098849 --db_config $HOME/repos/reportgen/tests/referral-db-config.json \
        #  --output /dev/stdout  --address_table_file reportgen/assets/addresses.csv
        return "compileMetadata" + \
               required('', self.blood_barcode) + \
               required('', self.tumor_barcode) + \
               required('--db_config ', self.referral_db_conf) + \
               required('--address_table_file ', self.addresses) + \
               required('--output ', self.output_json)


class CompileAlasccaGenomicJson(Job):
    def __init__(self, input_somatic_vcf, input_cn_calls, input_msisensor,
                 input_purity_qc, input_contam_qc, input_tcov_qc, input_ncov_qc,
                 output_json):
        Job.__init__(self)
        self.input_somatic_vcf = input_somatic_vcf
        self.input_cn_calls = input_cn_calls
        self.input_msisensor = input_msisensor
        self.input_purity_qc = input_purity_qc
        self.input_contam_qc = input_contam_qc
        self.input_tcov_qc = input_tcov_qc
        self.input_ncov_qc = input_ncov_qc
        self.output_json = output_json

    def command(self):
        return 'compileAlasccaGenomicReport ' + \
               required('', self.input_somatic_vcf) + \
               required('', self.input_cn_calls) + \
               required('', self.input_msisensor) + \
               required('--tumorCovJSON ', self.input_tcov_qc) + \
               required('--normalCovJSON ', self.input_ncov_qc) + \
               required('--purityJSON ', self.input_purity_qc) + \
               required('--contaminationJSON ', self.input_contam_qc) + \
               required('--output ', self.output_json)


class WriteAlasccaReport(Job):
    def __init__(self, input_genomic_json, input_metadata_json, output_pdf, alascca_only=False):
        Job.__init__(self)
        self.input_metadata_json = input_metadata_json
        self.input_genomic_json = input_genomic_json
        self.output_pdf = output_pdf
        self.alascca_only = alascca_only

    def command(self):
        tmpdir = "{}/write-alascca-report-{}".format(self.scratch, uuid.uuid4())
        mkdir_tmp_cmd = "mkdir -p {}".format(tmpdir)
        tmp_pdf = os.path.join(tmpdir, 'Report.pdf')
        cmd = 'writeAlasccaReport ' + \
              required(' --tmp_dir ', tmpdir) + \
              required(' --output_dir ', tmpdir) + \
              conditional(self.alascca_only, " --alascca_only ") + \
              required('', self.input_genomic_json) + \
              required('', self.input_metadata_json)

        cp_cmd = "cp {} {}".format(tmp_pdf, self.output_pdf)
        rmdir_cmd = "rm -r {}".format(tmpdir)

        return " && ".join([mkdir_tmp_cmd, cmd, cp_cmd, rmdir_cmd])
