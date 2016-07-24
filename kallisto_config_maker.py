import os
import sys
import json


class fastq_info:
    def __init__(self, directory):
        self.dir = directory
        self.all_fastqs = self.get_fastqs_list(self.dir)
        self.sample_fastq_dict = self.group_fastqs_by_sample(self.all_fastqs)
        self.paired_end = self.check_for_pair(self.sample_fastq_dict)

    def get_fastqs_list(self, dir):
        files_in_dir = os.listdir(dir)
        path = os.path.abspath(dir)
        fastq_files = [os.path.join(path, f) for f in files_in_dir
                       if 'fastq' in f]
        return fastq_files

    def group_fastqs_by_sample(self, fastq_files):
        sample_fastq_dict = {}
        for x in fastq_files:
            sample_name = x.split("/")[-1].split('_')[0]
            if sample_name in sample_fastq_dict:
                sample_fastq_dict[sample_name].append(x)
            else:
                sample_fastq_dict[sample_name] = [x]
        return sample_fastq_dict

    def check_for_pair(self, sample_fastq_dict):
        paired = []
        for sample in sample_fastq_dict:
            sample_pair = False
            for fastq in sample_fastq_dict[sample]:
                if 'R2' in fastq:
                    sample_pair = True
            paired.append(sample_pair)
        if all(paired):
            return True
        if not any(paired):
            return False
        raise AttributeError('Mix of paired and single end fastqs')


def check_for_gzip(sample_fastq_dict):
    for sample, fastq in sample_fastq_dict.items():
        if '.gz' in fastq:
            return True


def print_config_with_fastqs(fastq_class):
    if not os.path.exists('luigi.cfg'):
        with open('luigi.cfg', 'a') as f:
            print('[globalConfig]', file=f)
            print('exp_name:', file=f)
            print('#mnt_point:', file=f)
            print('fastq_dict:',
                  json.dumps(fastq_class.sample_fastq_dict), file=f)
            print('', file=f)
            print('[kallisto_index]', file=f)
            print('#transcript_fa:', file=f)
            print('[kallisto_quant]', file=f)
            print('paired:', fastq_class.paired_end, file=f)


def main():
    test = fastq_info(sys.argv[1])
    print_config_with_fastqs(test)

if __name__ == '__main__':
    main()
