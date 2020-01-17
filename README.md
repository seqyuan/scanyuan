# scanyuan
摘取了单细胞数据处理python包[scanpy](https://scanpy.readthedocs.io/en/stable/installation.html "scanpy")中的部分画图函数做了修改复现文章的图

# 安装
```
pip install scanyuan
```

# 使用

```python
import scanpy as sc
import scanyuan as scy

# 此处示例为读取loom文件，也可以是其他scanpy支持的数据格式
adata = sc.read_loom("/Users/yuanzan/Desktop/tmp/sdata.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
marker_genes = ['Stfa1', 'Ngp', 'Ccl5', 'Ccl4', 'BC100530', 'Gzma', 'Gata2', 'Cd74']

ax = scy.stacked_violin_t(adata, marker_genes, figsize=[8,7], groupby='ClusterName')
```

单细胞分析R包Seurat对象转loom格式请查看[僵小鱼的故事](http://www.seqyuan.com/jiangxiaoyu%20story.html)

# 效果
<img src="https://raw.githubusercontent.com/seqyuan/scanyuan/master/scy.png" width = "400" alt="" />


# 来自文献类的图
<img src="https://raw.githubusercontent.com/seqyuan/scanyuan/master/1111.png" width = "400" alt="" />

<img src="https://raw.githubusercontent.com/seqyuan/scanyuan/master/2222.png" width = "400" alt="" />



