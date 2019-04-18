# Pythonで作る2次元変位解析FEM
## 定式化
各要素の剛性マトリクス$Ke$は以下のようになる。
$$
K_{e}=\sum_{i=1}^4 \omega_iB^T(\xi_i,\eta_i)DB(\xi_i,\eta_i)det(J(\xi_i,\eta_i))
$$
また、ガウス・ルジャンドル求積により、
$$
K_{e}=\sum_{i=1}^4 \omega_iB^T(\xi_i,\eta_i)DB(\xi_i,\eta_i)det(J(\xi_i,\eta_i))
$$
また$B$マトリクスは以下のように表される。
$$
$$