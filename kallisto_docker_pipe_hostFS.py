import luigi
import subprocess
import os


class globalConfig(luigi.Config):
    exp_name = luigi.Parameter(default='test')
    mnt_point = luigi.Parameter(default='/home/jake/rna_data')
    fastq_dict = luigi.DictParameter(default={'test': 1})


class kallisto_index(luigi.Task):
    transcript_fa = luigi.Parameter(default='gencode.vM10.transcripts.fa.gz')

    def run(self):
        transcript_fa_path = '{data_store}/index/{transcript_fasta}'.format(
            data_store=globalConfig().mnt_point,
            transcript_fasta=self.transcript_fa
        )

        docker_cmd = [
            'docker', 'run',
            '-v', '{host_dir}:{mnt_point}'.format(
                host_dir=globalConfig().mnt_point,
                mnt_point=globalConfig().mnt_point),
            'kallisto', 'index',
            '-i', self.output().path,
            transcript_fa_path
            ]

        subprocess.call(docker_cmd)

    def output(self):
        index_name = self.transcript_fa.split(".fa")[0] + '.idx'
        return luigi.LocalTarget('{data_store}/index/{index_name}'.format(
            data_store=globalConfig().mnt_point,
            index_name=index_name
            ))


class kallisto_quant(luigi.Task):
    sample = luigi.Parameter()
    paired_end = luigi.BoolParameter(default=True)

    def requires(self):
        return kallisto_index()

    def run(self):
        try:
            os.makedirs('{data_store}/{exp_name}'.format(
                data_store=globalConfig().mnt_point,
                exp_name=globalConfig().exp_name)
            )
        except OSError:
            pass

        if self.paired_end:
            docker_cmd = [
                'docker', 'run',
                '-v', '{host_dir}:{mnt_point}'.format(
                    host_dir=globalConfig().mnt_point,
                    mnt_point=globalConfig().mnt_point),
                'kallisto', 'quant',
                '-i', self.input().path,
                '-o', self.output().path,
                ]
            docker_cmd += globalConfig().fastq_dict[self.sample]

        else:
            docker_cmd = [
                'docker', 'run',
                '-v', '{host_dir}:{mnt_point}'.format(
                    host_dir=globalConfig().mnt_point,
                    mnt_point=globalConfig().mnt_point),
                'kallisto', 'quant',
                '-i', self.input().path,
                '-o', self.output().path,
                '--single', '-l', '200', '-s', '20',
                ]
            docker_cmd += globalConfig().fastq_dict[self.sample]

        subprocess.call(docker_cmd)

    def output(self):
        return luigi.LocalTarget('{data_store}/{exp_name}/{sample}'.format(
            data_store=globalConfig().mnt_point,
            exp_name=globalConfig().exp_name),
            sample=self.sample)
