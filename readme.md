# Pythonで作る2次元変位解析FEM
## 定式化
各要素の剛性マトリクス$K_e$は以下のようになる。

$$
K_e=\iint_{\Omega_e} B^TDBdet(J)d\Omega_e
$$

また、ガウス・ルジャンドル求積により、

$$
K_e=\sum_{i=1}^4 \sum_{j=1}^4 \omega_i \omega_j B^T(\xi_i,\eta_i)DB(\xi_i,\eta_i)det(J(\xi_i,\eta_i))
$$
$$
\omega_i=\omega_j=1 \;\;(i,j=1,2,3,4)
$$
また$B$マトリクスは以下のように表される。

$$
B={B_1 }
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
B=AN_n=\left( \begin{array}{}
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