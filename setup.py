from setuptools import (
    setup,
)

setup(
    name="scanyuan",
    version = "0.0.5",
    author="seqyuan",
    author_email='seqyuan@gmail.com',
    url="https://github.com/seqyuan/scanyuan",
    download_url = "https://codeload.github.com/seqyuan/scanyuan/zip/master",
    description="modified scanpy functions to plot AnnData",
    long_description="""modified scanpy functions to plot AnnData""",
    license="MIT",
    packages=['scanyuan'],
    extras_require = {
        'scanpy' : [ 'scanpy']
    }
)


