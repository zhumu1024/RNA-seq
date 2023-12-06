# 8. 合并表达矩阵与标准化

# 8.1 合并

rm(list=ls())
# 删除当前工作空间中的所有对象，即清空所有已存在的变量、函数
# 创建一个名为 list 的变量
# ls() 函数是 R 语言中用于列出当前工作空间中对象的函数，包括变量、函数、数据框等
# 在R语言中，rm() 函数用于移除（删除）对象，可以是变量、函数、数据框等


setwd("~/project/rat/output/HTseq")
# setwd() 是 R 中的一个函数，用于设置工作目录
# 双引号在这里类似bash中的作用

# 得到文件样本编号
files <- list.files(".", "*.count")
# <- 是 R 中的赋值操作符，将右侧的值赋给左侧的变量
#  list.files() 是 R 语言中用于列出指定目录中文件和子目录的函数。它返回一个字符向量，其中包含目录中的文件和子目录的名称。
# . 表示当前工作目录
# "*.count" 是一个通配符模式，表示只列出文件名符合 *.count 模式的文件。

f_lists <- list()
# list() 函数用于创建列表
# f_lists <- list() 将这个空列表赋值给了变量 f_lists

for(i in files){
    prefix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)
    #gsub() 函数是用于替换字符串中的模式的函数  
    # (_\\w+)?\\.count 正则表达式   
    # \\ 在此表示转义
    # （ ）捕获组
    # _下划线
    # \\w 表示字母或者数字 + 多个  ？ 非贪婪模式

    # "" 空字符串
    #  perl=TRUE 表示perl兼容的正则表达式

    f_lists[[prefix]] = i
    # 字典，创建映射
}

id_list <- names(f_lists)
# 将字典中的key给了id_list
data <- list()
count <- 0

for(i in id_list){
  count <- count + 1
  a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id",i))
  # read.table() 函数用于从文件中读取数据框
  # sep="\t" 表示数据文件中使用制表符（\t）作为字段的分隔符
  # col.names = c("gene_id",i) 设置数据框的列名
  # 双中括号 [[]] 用于从列表中提取单个元素，而单中括号 [] 用于提取子集或多个元素。

  data[[count]] <- a
}

# 合并文件
data_merge <- data[[1]]
for(i in seq(2, length(id_list))){
# seq(2, length(id_list)) 创建了一个从2到length(id_list)的整数序列

    data_merge <- merge(data_merge, data[[i]],by="gene_id")
    # merge() 函数用于按照指定的列将两个数据框进行合并
    # data_merge, data[[i]] 两个要指定的列
    # by="gene_id" 按照这个列来合并
}

write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)
# 将 data_merge 这个数据框写入 名为 merge.csv的csv文件
# quote = FALSE 表示在写入 CSV 文件时不要包含字段值的引号
# row.names = FALSE 表示不在文件中包含行号


8.2 数据标准化

# 得到基因长度
library(GenomicFeatures)
# R编程语言中加载GenomicFeatures库的命令


# 构建Granges对象
txdb <- makeTxDbFromGFF("rn6.gff" )
# 使用 makeTxDbFromGFF 函数，该函数用于从 GFF（General Feature Format）文件创建一个 TranscriptDb 数据库对象（txdb）
# TranscriptDb 数据库对象是 GenomicFeatures 包中用于存储基因组注释信息的一种数据结构。这个对象包含了有关基因、转录本和外显子等基因组特征的详细信息

# 查找基因的外显子
exons_gene <- exonsBy(txdb, by = "gene")
# 使用 exonsBy 函数从中提取外显子信息。这里 by = "gene" 表示按基因进行分组，即提取每个基因的外显子信息。结果是一个列表

# 计算总长度
# reduce()、width()是Irange对象的方法
gene_len <- list()
for(i in names(exons_gene)){
# names(exons_gene) 返回 exons_gene 列表的基因名

    range_info = reduce(exons_gene[[i]])
    # 对当前基因的外显子信息进行合并，使用 reduce 函数。这样可以将相邻的外显子合并为一个或多个区域，减少冗余的信息

    width_info = width(range_info)
    #  计算合并后的区域的宽度，即每个区域的长度,返回向量
    sum_len    = sum(width_info)
    # 将当前基因的总长度存储到 gene_len 列表中，以基因名为索引

    gene_len[[i]] = sum_len
}

# 或者写为lapply的形式(快很多)
gene_len <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
# lapply 是一个函数，它对列表（或向量）中的每个元素应用指定的函数

data <- t(as.data.frame(gene_len))
# as.data.frame 函数用于将一个对象转换为数据框。在这里，gene_len 是一个列表，每个元素都是一个基因的总长度。通过 as.data.frame(gene_len)，将这个列表转换为一个数据框
# t 函数用于转置矩阵或数据框。在这里，对先前得到的数据框进行转置，将原本的行变为列，列变为行。

# 写入文件
write.table(data, file = "rn6_gene_len.tsv", row.names = TRUE, sep="\t", quote = FALSE, col.names = FALSE)


#计算 FPKM TPM

#!R
# =========== RPKM =============

gene_len_file <- "rn6_gene_len.tsv"
count_file <- "samples.count"

gene_len <- read.table(gene_len_file, header = FALSE, row.name = 1)
# gene_len 被赋值为一个数据框
# row.name = 1 表示将文件的第一列作为行名 

colnames(gene_len) <- c("length")
# 使用 colnames 函数给数据框的列（唯一的列）命名

count <- read.table(count_file, header = FALSE, row.name = 1)
colnames(count) <- c("count")
# all read number
all_count <- sum(count["count"])
# count["count"]： 这部分代码使用方括号 [] 进行子集选择    选择count这一列
RPKM <- c()
# 创建一个空的向量 RPKM，用于存储计算得到的RPKM值

for(i in row.names(count)){
    count_ = 0
    exon_kb = 1
    rpkm = 0
    count_ = count[i, ] # 选取第i行
    # 从 count 数据框中选择第 i 行，即当前基因的reads数
    exon_kb  = gene_len[i, ] / 1000
    rpkm    = (10 ^ 6 * count_ ) / (exon_kb * all_count )
    RPKM = c(RPKM, rpkm)
}


# =========== 计算TPM ============
# 首先得到总的结果
sum_ <- 0
for(i in row.names(count)){
    count_ = 0
    exon = 1
    count_ = count[i, ]
    exon  = gene_len[i, ]
    value = count_ / exon
    if(is.na(value)){  # 检查 RPK 值是否为 NA（可能是由于分母为零等原因）
        print(paste(i, " is error! please check"))
    }else{
        sum_ = sum_ + value
    }
}

TPM <- c()
for(i in row.names(count)){
    count_ = 0
    exon = 1
    count_ = count[i, ]
    exon  = gene_len[i, ]
    tpm = (10 ^ 6 * count_ / exon ) / sum_
    TPM = c(TPM, tpm)
}

count["RPKM"] <- RPKM # count["RPKM"] <- RPKM： 将名为 "RPKM" 的列添加到数据框 count 中，并将 RPKM 向量的值赋给这一列。
count["TPM"] <- TPM       
           
write.table(count, "123.normalize.count", col.names = TRUE, row.names = TRUE, sep="\t", quote = FALSE)


# 9 差异表达分析

cd ~/project/rat/output/HTseq

cat merge.csv | grep -E "ENSRNOG00000018630|ENSRNOG00000034254"  ## E 扩展的正则表达式  这里面一般是多选

9.1 数据前处理

##删除HTseq-count的总结行

dataframe <- read.csv("merge.csv", header=TRUE, row.names = 1) 

# 去除前面5行
countdata <- dataframe[-(1:5),]  #使用负索引 - (1:5) 表示移除数据框的第1到第5行  通过负索引表示要排除的行

# 查看数据
head(countdata)  

# 得到行的名
row_names <- row.names(countdata) 

# 开始替换
name_replace <- gsub("\\.\\w+","", row.names(countdata))  # gsub  替换  \\. 因为点在正则表达式中有特殊意义，需要使用两个反斜杠 \\ 来表示实际的点

row.names(countdata) <- name_replace

#去除低表达的基因
countdata <- countdata[rowSums(countdata) > 0,]   # rowSums 函数计算数据框 countdata 中每行的元素之和  
# countdata[rowSums(countdata) > 0,]： 使用逻辑索引   保留大于0的行


9.2  差异表达分析
# 使用bioconductor进行安装
source("http://bioconductor.org/biocLite.R") #  source() 函数： 用于在 R 中执行脚本文件
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/") # options() 函数： 用于设置和查询 R 的全局选项  这里是使用镜像文件

# 安装包
biocLite("DESeq2")  # 用于差异表达分析的工具
biocLite("pheatmap") # 一个用于创建热图的包
biocLite("biomaRt") # 提供了访问 BioMart 数据库的接口
biocLite("org.Rn.eg.db") # 这是一个包含了Rattus norvegicus（大鼠）物种的基因注释信息的数据库
biocLite("clusterProfiler") # 这是用于富集分析和聚类分析的工具

# 加载
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(org.Rn.eg.db)
library(clusterProfiler)


9.2.2 构建对象

cat <<EOF >./phenotype/phenotype.csv  #Here Document（
"ids","state","condition","treatment"
"SRR2240185","Liver cirrhosis","DEN","treatment"
"SRR2240186","Liver cirrhosis","DEN","treatment"
"SRR2240187","Healthy control","PBS","control"
"SRR2240228","Healthy control","PBS","control"
EOF

# 刚才countdata已经得到
countdata

# 读取样本分组信息(注意，需要加上row.names = 1, header = TRUE，将行列名需要看好)
coldata <- read.table("../phenotype/phenotype.csv", row.names = 1, header = TRUE, sep = "," ) 
# 确认一下行列名是否有（不是简单的数值）
head(coldata)
# 调整数据顺序
countdata <- countdata[row.names(coldata)]  # countdata[] 表示对数据框（或矩阵）countdata 进行子集选择（subset）的操作  

# 构建dds对象
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~ treatment)  
# 使用 DESeq2 包中的函数 DESeqDataSetFromMatrix 创建一个差异表达分析的数据集对象（DESeqDataSet）
# 参数 countData 指定了包含基因计数数据的矩阵或数据框
# 参数 colData 指定了包含样本信息的数据框，其中包括样本的条件、处理等信息
# 参数 design 指定了线性模型的设计。在这里，使用了 "~ treatment"，表示在模型中考虑处理（treatment）作为差异表达的主要因素。


# 查看dds
dds

9.2.3 样本相关性

# 接续着上面的构建得到的dds对象
# DEseq2包提供了相应的函数
vsdata <- rlog(dds, blind=FALSE) # rlog 函数用于进行对数变换，以减小基因表达数据中的离散性和异方差性。 
# blind=FALSE 参数表示在进行对数变换时，不采用“盲变换”策略。盲变换是一种在不使用先验信息的情况下，
# 通过对数变换来调整数据的方法。将 blind 设置为 FALSE 意味着考虑到样本之间的差异性。
# intgroup 分组
plotPCA(vsdata, intgroup="treatment") + ylim(-10, 10)
# plotPCA 是 DESeq2 包中的一个函数，用于绘制主成分分析（PCA）的图形。
# 在基因表达数据分析中，PCA常用于可视化样本之间的整体差异，帮助我们观察样本是否存在聚类趋势，以及哪些主成分方向对于解释数据的变异贡献较大。
# vsdata 是之前通过 rlog 转换后的数据，intgroup="treatment" 指定了用于着色数据点的分组信息

#+ ylim(-10, 10)： 这部分对绘制的PCA图进行修饰。ylim(-10, 10) 设置了y轴的显示范围，限制了y轴的取值范围在-10到10之间




#sample-to-sample distances热图
# 颜色管理包（不是必须）
library("RColorBrewer")

# 得到数据对象中基因的计数的转化值
gene_data_transform <- assay(vsdata) #assay 函数主要用于从 SummarizedExperiment 或 DESeqDataSet 对象中提取矩阵型数据。
# 使用t()进行转置
# 使用dist方法求样本之间的距离--衡量样本之间的相似性或差异性
sampleDists <- dist(t(gene_data_transform))

# 转化为矩阵用于后续pheatmap()方法的输入
sampleDistMatrix <- as.matrix(sampleDists) # as.matrix 函数用于将一个R对象转换为矩阵。这个对象可以是数据框、数组、列表等。

# 将矩阵的名称进行修改
# rownames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")
# colnames(sampleDistMatrix) <- paste(vsdata$treatment, vsdata$condition, vsdata$ids, sep="-")

# 设置色盘
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# brewer.pal 是 RColorBrewer 包中的一个函数，用于获取预定义的颜色调色板
# rev 函数是用于反转向量的函数。在这里，它被用于反转由 brewer.pal 生成的颜色向量，以便颜色从浅到深排列
# colorRampPalette 是一个函数，用于创建颜色渐变的色盘。它接受一个颜色向量，并返回一个能够生成渐变颜色的函数。
# (255) 表示调用先前生成的颜色渐变色盘函数，以获得包含 255 种颜色的向量


# 绘制热图与聚类
pheatmap(sampleDistMatrix, # 一个矩阵，表示样本之间的距离。通常是通过欧氏距离或其他距离度量计算得到的距离矩阵。
         clustering_distance_rows=sampleDists, # 用于指定在行（样本）上应用的聚类方法的距离矩阵。在这里，使用了之前计算得到的 sampleDists
         clustering_distance_cols=sampleDists, # 用于指定在列上应用的聚类方法的距离矩阵。同样使用了之前计算得到的 sampleDists
         col=colors)


# 9.2.4 差异基因

#使用DESeq()方法计算不同组别间的基因的表达差异，它的输入是上一步构建的dss数据对象
# 改变样本组别顺序
dds$treatment <- factor(as.vector(dds$treatment), levels = c("control","treatment"))
# #目的是重新排序 dds 中的 treatment 列，确保"control"在前，"treatment"在后
# as.vector 用于将 R 中的数据结构转换为向量。在 R 中，数据结构如列表、矩阵等有时候不是向量，使用 as.vector 可以将其强制转换为向量。
# factor 用于创建因子型变量  

# 基于统计学方法进行计算
dds <- DESeq(dds)
#DESeq 函数的主要目标是使用负二项分布模型估计基因的表达差异，并进行假设检验，计算 p 值等。执行这个函数后，dds 对象中将包含有关每个基因的差异表达分析的结果

#查看实验组与对照组的对比结果
result <- results(dds, pAdjustMethod = "fdr", alpha = 0.05) # FDR（False Discovery Rate，假发现率）进行多重检验校正  # 显著性水平取0.05
#results 函数用于从已经进行差异表达分析的 DESeqDataSet 对象中提取结果

# 查看结果
head(result)

#将结果按照p-value进行排序
result_order <- result[order(result$pvalue),] #  R 中的排序函数，用于按照某个变量的值对结果进行排序   $pvalue 按照这个来排序
head(result_order)

#总结基因上下调情况
summary(result_order)  #summary 函数用于生成一个对结果的简要摘要 如上调基因和下调基因的数量等

# 查看显著的基因数量
table(result_order$padj<0.05)  # table 函数用于统计逻辑向量中的各个值的频数。在这里，result_order$padj < 0.05 是一个逻辑向量，

# 将数据保存起来

# 新建文件夹
dir.create("../DESeq2")
# 不用按照padj排序的结果，就保存按照基因名排序的
write.csv(result, file="../DESeq2/results.csv", quote = F)


----------------------------------------------------------------------------------------------------------------------------------------
##10. 提取差异表达基因与注释

# padj 小于 0.05 并且 Log2FC 大于 1 或者小于 -1
diff_gene <- subset(result_order, padj < 0.05 & abs(log2FoldChange) > 1)  逻辑部分

# 查看数据框的大小
dim(diff_gene) # dim() 是 R 语言中的一个用于获取对象维度的函数。对于矩阵或数据框等二维结构，dim() 函数返回一个包含行数和列数的向量。

# 把差异基因写入到文件中
dir.create("../DESeq2/")
write.csv(diff_gene, file="../DESeq2/difference.csv", quote = F) 

#10.2 使用ClusterProfiler对基因的ID进行转化

# 首先安装ClusterProfiler
source("http://bioconductor.org/biocLite.R")
# 安装clusterProfiler包
biocLite("clusterProfiler")
# 一个 R 语言的生物信息学包，用于进行生物学数据的富集分析和聚类可视化
# 
# 这里我们分析的是大鼠，安装大鼠的数据库
biocLite("org.Rn.eg.db")

# 加载包
library(clusterProfiler)
library(org.Rn.eg.db)

# 得到基因ID(这个ID是Ensembl数据库的编号)
ensembl_gene_id <- row.names(diff_gene)

# 转换函数
ensembl_id_transform <- function(ENSEMBL_ID){
    # geneID是输入的基因ID，fromType是输入的ID类型，toType是输出的ID类型，OrgDb注释的db文件，drop表示是否剔除NA数据
    a = bitr(ENSEMBL_ID, fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Rn.eg.db")
    #，bitr 函数的作用是将 Ensembl 基因 ID 转换为基因符号和 Entrez 基因 ID
    # "SYMBOL"： 表示希望将输入的基因 ID 转换为基因符号（Gene Symbol）
    # "ENTREZID"： 表示希望将输入的基因 ID 转换为 Entrez 基因 ID   ncbi
    return(a)
    #a 是一个包含了从 Ensembl 基因 ID 转换得到的基因符号（Symbol）和 Entrez 基因 ID 的数据框。
    # return 用于指定函数的返回值
}

# 开始转化
ensembl_id_transform(ensembl_gene_id)

# 10.3 使用biomaRt进行注释
# 安装
biocLite("biomaRt")
#biomaRt 是一个 R 语言的生物信息学包，用于访问和查询 BioMart 数据库。BioMart 是一个集成了多个生物学数据库的数据查询工具，允许用户从不同数据库中检索和下载生物学数据

# 加载
library("biomaRt")

# 选择数据库
mart <- useDataset("rnorvegicus_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
#useDataset 函数用于指定要使用 BioMart 数据库中的哪个数据集
# useMart 函数用于指定使用哪个 BioMart 数据库
#"ENSEMBL_MART_ENSEMBL" 表示使用 Ensembl 数据库。这是一个 Ensembl 提供的 BioMart 数据库，包含了有关多个物种的基因、转录本和蛋白质的注释信息。
# "rnorvegicus_gene_ensembl" 表示选择大鼠（Rattus norvegicus）的基因注释数据集

# 得到基因ID(这个ID是Ensembl数据库的编号)
ensembl_gene_id <- row.names(diff_gene) 
rat_symbols <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id", "description"), 
#getBM 函数，通过 mart 对象从 BioMart 数据库中检索大鼠基因的注释信息
#getBM 函数用于从 BioMart 数据库中获取特定属性（attributes）的注释信息  attributes 参数指定你想要检索的属性
# 
filters = 'ensembl_gene_id', values = ensembl_gene_id, mart = mart)
#  filters 参数指定过滤器的类型
# values 参数指定过滤器的取值，即之前提取的 Ensembl 基因 ID
# mart 参数指定了使用哪个 mart 对象，即之前创建的 BioMart 数据库对象

#将基因矩阵与symbols合并
# 生成用于合并的列
diff_gene$ensembl_gene_id <- ensembl_gene_id
# 将DESeq2对象转换为数据库
diff_gene_dataframe <- as.data.frame(diff_gene) #as.data.frame 是一个 R 语言中的函数，用于将其他数据结构转换为数据框（data frame）。
# 合并
diff_gene_symbols <- merge(diff_gene_dataframe, rat_symbols, by = c("ensembl_gene_id"))
# merge 函数，将两个数据框 diff_gene_dataframe 和 rat_symbols 按照共同的列 "ensembl_gene_id" 进行合并

# 将数据存储起来
write.table(result, "../stat/all_gene.tsv", sep="\t", quote = FALSE)
write.table(diff_gene_symbols, "../stat/diff_gene.tsv", row.names = F,sep="\t", quote = FALSE)

# 统计样本的差异基因
echo -e "sample\tnum" > all_samples.tsv
for i in $(ls);
do
    if [ -d ${i} ];
    then
        prefix=$i
        diff_num=$(cat $i/diff_gene.tsv | tail -n+2 | wc -l)
        #tail -n+2 跳过文件的第一行，即表头   wc -l 统计剩余的行数。
        echo -e "${prefix}\t${diff_num}" >> all_samples.tsv
        # 
    fi
done

# 使用R绘图
library(ggplot2)
data <- read.table("all_samples.tsv", header = T)

pdf("samples_diff_gene_num.pdf") #这行代码打开一个 PDF 设备，将后续绘制的图形保存为 PDF 文件
  ggplot(data=data, aes(x=sample, y=num, fill=sample)) + #data=data 指定了用于绘图的数据框  
  #aes(x=sample, y=num, fill=sample) 是一个映射，指定了横坐标、纵坐标和颜色填充的映射关系。
  geom_bar(stat = "identity", position = "dodge") + #ggplot2 中的几何对象，用于添加柱状图 
  # stat = "identity" 表示使用数据中的实际值作为柱状图的高度
  # position = "dodge" 表示将多个柱状图分开显示，避免它们重叠在一起
  labs(x = "samples",y = "num",title = "different gene number")
  # ggplot2 中的函数，用于设置图形的标签
  # x = "samples" 设置横坐标的标签为 "samples"。
  # y = "num" 设置纵坐标的标签为 "num"
dev.off()





