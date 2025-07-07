# weyl
第一原理計算と対称性計算の結合を目指す
- spglib.spglib.get_symmetry()
- spglib.spglib.SpglibDataset
- spglib.spglib.get_symmetry_dataset()
## 手順
- 対称操作で波動関数の指標表入手
- 対称性をベクトルで表す
- 可能なら、
  - 点群の既約表現をLLMに生成させる
  - 点群の操作とspglibの操作を対応づけて既約表現の記号を得る
