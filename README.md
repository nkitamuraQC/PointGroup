# weyl
第一原理計算と対称性計算の結合を目指す
- spglib.spglib.get_symmetry()
- spglib.spglib.SpglibDataset
- spglib.spglib.get_symmetry_dataset()
## 手順
- 対称操作で波動関数の指標表入手
- 対称性を指標ベクトルで表す(<-これで十分)
- SALCを作る
- 可能なら、
  - 点群の既約表現をLLMに生成させる
  - 点群の操作とspglibの操作を対応づけて既約表現の記号を得る
## 群論メモ
- 正準分子軌道は、群の既約表現に従うように構築されている
- 異なる既約表現の指標の内積はゼロ
