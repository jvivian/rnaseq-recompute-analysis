from setuptools import setup

setup(name='rnaseq_recompute_analysis',
      version='1.0a1',
      description='Repository for reproducing analyses for the CGL RNA-seq recompute',
      url='http://github.com/jvivian/rnaseq-recompute-analysis',
      author='John Vivian',
      author_email='jtvivian@gmail.com',
      license='MIT',
      packages=['recompute_analysis'],
      install_requires=['tqdm', 'pandas', 'numpy'],
      entry_points = {
          'console_scripts': [
              'create-project = recompute_analysis.create_project:main',
              'recompute-analysis = recompute_analysis.recompute_analysis:main'
          ]
      })
