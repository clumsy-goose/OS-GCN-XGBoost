# 基于机器学习的水稻抗旱基因研究

## dataSource:原始数据

## dataPreprocessing:数据预处理

* analyze:数据分析

* processing:数据处理
* result:中间结果数据

## dataProcessedResult:数据预处理的结果

* geneWeight:由基因互作网络提取的基因权重
* phenotype_rpkm:性状表达-基因表达关系

## model:模型

* GCN: GCN处理基因互作网络，获得基因权重特征
* XGBoost：使用XGBost拟合性状表达-基因表达关系
* GCN+XGBoost：基因权重特征结合性状表达-基因表达关系，使用XGBoost拟合
