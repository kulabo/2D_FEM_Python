# Pythonで作る2次元変位解析FEM
## 定式化
弾性体の支配方程式は、いわゆるフックの法則である以下の式で表される。
$$
F=KU
$$

この式は解析対象全体についての式だが、これを有限要素に分解(離散化)する。1要素についての支配方程式は以下である。
$$
F_e=K_eU_e
$$

ただし、各要素の剛性マトリクス$K_e$は
$$
K=\sum_e^{AllElements}K_e
$$

であり、$K_e$は以下のようになる。
$$
K_e=\iint_{\Omega_e} B^TDBdet(J)d\Omega_e
$$

また、ガウス・ルジャンドル求積により、
$$
K_e=\sum_{i=1}^4 \sum_{j=1}^4 \omega_i \omega_j B^T(\xi_{i,j},\eta_{i,j})DB(\xi_{i,j},\eta_{i,j})det(J(\xi_{i,j},\eta_{i,j}))
$$

$$
\omega_i=\omega_j=1 \;\;(i,j=1,2,3,4)
$$

$$
(\xi_{i,j},\eta_{i,j})=
(\pm \frac{1}{\sqrt 3},\pm \frac{1}{\sqrt 3}),
(\pm \frac{1}{\sqrt 3},\mp \frac{1}{\sqrt 3})
$$

また$B$マトリクスは以下のように表される。

$$
B=\left( \begin{array}{}
B_1 & B_2 & B_3 & B_4 
\end{array} \right)
$$

$$
\epsilon=\left( \begin{array}{} 
\epsilon_x \\
\epsilon_y \\
\gamma_{xy}
\end{array} \right)
=\left( \begin{array}{}
\frac{\partial u}{\partial x} \\
\frac{\partial v}{\partial y} \\
\frac{\partial v}{\partial x}+
\frac{\partial u}{\partial y}
\end{array} \right)
=\left( \begin{array}{}
\frac{\partial}{\partial y} & 0 \\
0 & \frac{\partial}{\partial x} \\
\frac{\partial}{\partial x} & \frac{\partial}{\partial y}
\end{array} \right)
\left( \begin{array}{}
u \\
v
\end{array} \right)
=Au
$$

$$
D=\frac{E}{1-\mu^2}
\left( \begin{array}{}
1 & \mu & 0\\ 
\mu & 1 & 0 \\
0 & 0 & \frac{1-\mu}{2} 
\end{array} \right)
$$

$$
B_n=AN_n=\left( \begin{array}{}
\frac{\partial}{\partial y} & 0 \\
0 & \frac{\partial}{\partial x} \\
\frac{\partial}{\partial x} & \frac{\partial}{\partial y}
\end{array} \right)
\left( \begin{array}{}
N_n & 0 \\ 
0 & N_n
\end{array} \right)
=\left( \begin{array}{}
\frac{\partial N_n}{\partial y} & 0 \\
0 & \frac{\partial N_n}{\partial x} \\
\frac{\partial N_n}{\partial x} & \frac{\partial N_n}{\partial y}
\end{array} \right)
$$

$$
\frac{\partial N_n}{\partial \xi}=
\frac{\partial x}{\partial \xi} \frac{\partial N_n}{\partial x}+
\frac{\partial y}{\partial \xi} \frac{\partial N_n}{\partial y}
$$

$$
\frac{\partial N_n}{\partial \eta}=
\frac{\partial x}{\partial \eta} \frac{\partial N_n}{\partial x}+
\frac{\partial y}{\partial \eta} \frac{\partial N_n}{\partial y}
$$

より、

$$
\left( \begin{array}{}
\frac{\partial N_n}{\partial \xi} \\ 
\frac{\partial N_n}{\partial \eta} 
\end{array} \right)=
\left( \begin{array}{}
\frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\ 
\frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta}
\end{array} \right)
\left( \begin{array}{}
\frac{\partial N_n}{\partial \xi} \\ 
\frac{\partial N_n}{\partial \eta} 
\end{array} \right)=
J
\left( \begin{array}{}
\frac{\partial N_n}{\partial \xi} \\ 
\frac{\partial N_n}{\partial \eta} 
\end{array} \right)
$$

$$
x = \sum_{n=1}^4 N_n x_n ,\;\;
y = \sum_{n=1}^4 N_n y_n
$$

$$
\frac{\partial x}{\partial \xi} = \sum_{n=1}^4 \frac{\partial N_n}{\partial \xi} x_n
,\;\;
\frac{\partial x}{\partial \eta} = \sum_{n=1}^4 \frac{\partial N_n}{\partial \eta} x_n
$$

$$
\frac{\partial y}{\partial \xi} = \sum_{n=1}^4 \frac{\partial N_n}{\partial \xi} y_n
,\;\;
\frac{\partial y}{\partial \eta} = \sum_{n=1}^4 \frac{\partial N_n}{\partial \eta} y_n
$$

$$
J=
\left( \begin{array}{}
\frac{\partial N_1}{\partial \xi} & 
\frac{\partial N_2}{\partial \xi} &
\frac{\partial N_3}{\partial \xi} &
\frac{\partial N_4}{\partial \xi} \\ 
\frac{\partial N_1}{\partial \eta} &
\frac{\partial N_2}{\partial \eta} &
\frac{\partial N_3}{\partial \eta} &
\frac{\partial N_4}{\partial \eta} 
\end{array} \right)
\left( \begin{array}{}
x_1 & y_1 \\ 
x_2 & y_2 \\ 
x_3 & y_3 \\ 
x_4 & y_4 
\end{array} \right)
$$

$$
\left( \begin{array}{}
N_1 \\ 
N_2 \\ 
N_3 \\ 
N_4 
\end{array} \right)=
\frac{1}{4}
\left( \begin{array}{}
(1-\xi)(1-\eta) \\ 
(1+\xi)(1-\eta) \\
(1+\xi)(1+\eta) \\
(1-\xi)(1+\eta) \\
\end{array} \right)
$$

$$
\frac{\partial N_1}{\partial \xi}=-\frac{1-\eta}{4} ,\;\;
\frac{\partial N_2}{\partial \xi}=\frac{1-\eta}{4} ,\;\;
\frac{\partial N_3}{\partial \xi}=\frac{1+\eta}{4} ,\;\;
\frac{\partial N_4}{\partial \xi}=-\frac{1+\eta}{4}
$$

$$
\frac{\partial N_1}{\partial \eta}=-\frac{1-\xi}{4} ,\;\;
\frac{\partial N_2}{\partial \eta}=-\frac{1+\xi}{4} ,\;\;
\frac{\partial N_3}{\partial \eta}=\frac{1+\xi}{4} ,\;\;
\frac{\partial N_4}{\partial \eta}=\frac{1-\xi}{4}
$$