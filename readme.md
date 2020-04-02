# Pythonで作る2次元変位解析FEM

Pythonを用いた2次元弾性体に対する変位解析のデモプログラムです。コードの解説については[こちらのサイト](https://kulabo.github.io/2D_FEM_Python/python-2d-fem-01.html)をご確認ください。


## 1.はじめに
### 1.1 Pythonとは
Pythonは、いわゆるLL(Lightweight Language)のひとつである。ここでのLightweightというのは、動作が軽量という意味ではなく、習得・学習・使用が容易(軽量)であるということから付けられている。Pythonはその学習のしやすさ、それに伴うユーザーの多さ、パッケージの豊富さから、近年のAIブームにおける使用言語のデファクトスタンダードとなっている。もしあなたが今Pythonを使ったことがなかったとしても、これを機にPythonを始めてみることをおすすめする。

### 1.2 有限要素法とは
有限要素法(FEM : Fenite Element Method)とは、解析したい物や場を細かい要素に区切り、その細かい要素について計算した結果を重ね合わせて全体を求める手法である。流れ場の解析や構造物の応力解析などに広く使われる手法である。本レポートでは、2次元の構造物を対象とした応力解析問題について有限要素法を用いる。

## 2.定式化
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

一般に積分をプログラム上で行うことは困難であるから、数値積分によって計算し易い形式に近似する。数値積分法の1つであるガウス・ルジャンドル求積を適用して、
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

となる。$(\xi_{i,j},\eta_{i,j})$はガウス・ルジャンドル求積の積分点の正規化座標での座標値である。
今、上式においての未知数は$B$、$D$、$J$である。それらの具体形を以下に示す。
ひずみ$\epsilon$は以下のように表される。
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
また$B$マトリクスは上記の$A$を用いて以下のように表される。

$$
B=\left( \begin{array}{}
B_1 & B_2 & B_3 & B_4 
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
(n=1,2,3,4)
$$

ここで、$B_n$の成分は偏微分の連鎖則から以下のように表される。
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

以上の関係を行列形式で表せば、

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

となる。ここで$J$はヤコビマトリクスである。

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

<img src="https://latex.codecogs.com/gif.latex?\ \frac{\partial N_1}{\partial \xi}=-\frac{1-\eta}{4} ,\;\;\frac{\partial N_2}{\partial \xi}=\frac{1-\eta}{4} ,\;\;" />

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

$D$マトリクスは以下である。これは平面応力状態を仮定した際の式である。平面応力状態は薄い平板が面方向に荷重を受ける状態であり、板に垂直方向($z$方向)の応力を0とみなす仮定である。

$$
D=\frac{E}{1-\nu^2}
\left( \begin{array}{}
1 & \nu & 0\\ 
\nu & 1 & 0 \\
0 & 0 & \frac{1-\nu}{2} 
\end{array} \right)
$$

$E$はヤング率、$\nu$はポアソン比である。
参考として、平面ひずみ状態での$D$マトリクスを以下に示す。平面ひずみ状態では、$z$方向を非常に長いと仮定して**$z$方向のひずみ**を0とみなす。

$$
D=\frac{E}{(1+\nu)(1-2\nu)}
\left( \begin{array}{}
1-\nu & \nu & 0\\ 
\nu & 1-\nu & 0 \\
0 & 0 & \frac{1-2\nu}{2} 
\end{array} \right)
$$

## 3.実装
### 概要
本プログラムはトポロジー最適化という設計最適化手法のためのFEMとして設計した関係から、拘束条件や荷重条件、解析領域が変化せず、各要素のヤング率だけが変化するような繰り返しの解析がしやすくなっている。
FEMをクラスとして実装し、コンストラクタに荷重条件などの変化しない値を渡し、そのインスタンスの持つ計算メソッドに変化するヤング率(実際には無次元密度)を渡してループさせて最適化を行なう。

まずは大まかにFEMクラス