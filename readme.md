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
x = \sum_{n=1}^4 N_n x_n
$$
$$
y = \sum_{n=1}^4 N_n y_n
$$
$$
y = \sum_{n=1}^4 \frac{\partial N_n}{\partial \xi} x_n
$$