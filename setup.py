from setuptools import setup, find_packages

setup(name='rnaseq_recompute_analysis',
      version='1.0a1',
      description='Repository for reproducing analyses for the CGL RNA-seq recompute',
      url='http://github.com/jvivian/rnaseq-recompute-analysis',
      author='John Vivian',
      author_email='jtvivian@gmail.com',
      license='MIT',
      packages=find_packages(),#['recompute_analysis'],
      install_requires=['tqdm',
                        'pandas',
                        'numpy',
                        'synapseclient',
                        'matplotlib',
                        'sklearn',
                        'scipy',
                        'futures',
                        'toil>=3.6.0'],
      entry_points = {
          'console_scripts': [
              'create-project = create_project:main',
              'recompute-analysis = recompute_analysis:main',
              'pairwise-deseq2 = toil_pipelines.pairwise_deseq2:main'
          ]
      })
