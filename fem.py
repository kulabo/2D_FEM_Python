import time
import numpy as np
import matplotlib.pyplot as plt


class FEM:
    def __init__(self, pnl, nx, ny, mesh_size, fix_nodes, F):
        # self.rho = rho
        self.pnl = pnl
        self.nx = nx
        self.ny = ny
        self.mesh_size = mesh_size
        self.fix_nodes = fix_nodes
        self.F = F

        #self.E0 = 2.1*(10**11)
        self.E0 = 1.0
        self.Emin = 1e-3

        self.all_element_count = nx*ny
        self.all_node_count = (nx+1)*(ny+1)
        self.all_vector_count = 2*(nx+1)*(ny+1)

        self.free_nodes = np.setdiff1d(
            np.arange(self.all_vector_count), fix_nodes)

        self.node_coordinate_values = np.array(
            [[x, y] for y in range(ny+1) for x in range(nx+1)])*mesh_size

        # start from 1
        # [[elem1node1,elem1node2,elem1node3,elem1node4],
        #  [elem2node1,elem2node2,elem2node3,elem2node4],...]
        '''
        node numbers for each elements
        4 --- 3
        |     |
        |     |
        1 --- 2
        '''
        self.node_connection = np.array(
            [[(i+(nx+1)*j), (i+(nx+1)*j)+1, (i+(nx+1)*j)+nx+2, (i+(nx+1)*j)+nx+1] for i in range(1, nx+1) for j in range(ny)])

        self.point = np.array([[-3 ** (-0.5), -3 ** (-0.5)],
                               [3 ** (-0.5),  -3 ** (-0.5)],
                               [3 ** (-0.5),  3 ** (-0.5)],
                               [-3 ** (-0.5),  3 ** (-0.5)]])
        self.weight = np.array([1, 1, 1, 1])

        nu = 0.3
        # nu: poisson ratio
        # 平面応力状態
        self._Cmat = np.array([[1, nu, 0],
                               [nu, 1, 0],
                               [0, 0, (1-nu)/2]]) / (1-nu ** 2)
        # Dmat=E*Cmat

    # fenite element method
    def fem(self, rho):
        self.rho = rho

        K = self._Kmat()
        # np.savetxt('K.csv', K, delimiter=',')

        K_free = np.stack([K[self.free_nodes, col_index]
                           for col_index in self.free_nodes]).T
        U = np.zeros(self.all_vector_count)
        U[self.free_nodes] = np.linalg.solve(K_free, self.F[self.free_nodes])

        self.U = U

        # l: mean compliance = UtKU
        l = (U.T @ K) @ U
        return U, l

    # stiffness matrix
    def _Kmat(self):
        K = np.zeros([self.all_vector_count, self.all_vector_count])
        for y in range(self.ny):
            for x in range(self.nx):
                Ke = self._Kemat(x, y)
                top1 = 2*((self.ny+1)*x+y)
                top2 = 2*((self.ny+1)*(x+1)+y)
                K[top1:top1+4, top1:top1+4] += Ke[0:4, 0:4]
                K[top1:top1+4, top2:top2+4] += Ke[0:4, 4:8]
                K[top2:top2+4, top1:top1+4] += Ke[4:8, 0:4]
                K[top2:top2+4, top2:top2+4] += Ke[4:8, 4:8]
        
        #print(Ke==K)
        return K

    # stiffness matrix of one element
    def _Kemat(self, x, y):
        elem_rho = self.rho[x, y]
        Ke = np.zeros((8, 8))

        E = (self.E0-self.Emin)*elem_rho ** self.pnl+self.Emin
        D = E*self._Cmat

        for n in range(4):
            xi, eta = self.point[n]
            w = self.weight[n]
            dNdxi = self._dNdxi(eta)
            dNdeta = self._dNdeta(xi)
            J = self._Jmat(x, y, dNdxi, dNdeta)
            B = self._Bmat(J, dNdxi, dNdeta)
            # B = self._Bmat(dNdxi, dNdeta)
            print(B)
            # Bマトリクスは算出方法が違うので値は同じにならない
            Ke += w*(B.T @ D) @ B*np.linalg.det(J)
            # print(Ke)

        return Ke

    def _Jmat(self, x, y, dNdxi, dNdeta):
        '''
        node numbers for each elements
        4 --- 3
        |     |
        |     |
        1 --- 2
        '''

        x_nodes = np.array([x, x+1, x+1, x])
        y_nodes = np.array([y, y, y+1, y+1])

        J = np.stack([dNdxi, dNdeta]) @ np.stack([x_nodes, y_nodes]).T
        return J

    # partial derivative of shape function (N) with respect to xi
    def _dNdxi(self, eta):
        return np.array([-0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta)])

    # partial derivative of shape function (N) with respect to eta
    def _dNdeta(self, xi):
        return np.array([-0.25*(1-xi), -0.25*(1+xi), 0.25*(1+xi), 0.25*(1-xi)])

    
    def _Bmat(self, J, dNdxi, dNdeta):
        B = np.zeros([3, 8])
        dNdx_dNdy = np.linalg.inv(J) @ np.stack([dNdxi, dNdeta])
        for i in range(4):
            dNidx = dNdx_dNdy[0][i]
            dNidy = dNdx_dNdy[1][i]

            B[0][2*i] = dNidx
            B[1][2*i+1] = dNidy
            B[2][2*i] = dNidy
            B[2][2*i+1] = dNidx
        return B
    

    def plot_mesh(self):
        x_list = [row[0] for row in self.node_coordinate_values]
        y_list = [row[1] for row in self.node_coordinate_values]
        #x_list[2], x_list[3] = x_list[3], x_list[2]
        #y_list[2], y_list[3] = y_list[3], y_list[2]

        #plt.fill(x_list, y_list, c="r",alpha=0.5)
        print(x_list,y_list)
        plt.scatter(x_list, y_list, c="r", alpha=0.5)
        #plt.plot(x_list, y_list)

        '''
        x_list[0] += self.U[2*0]
        x_list[1] += self.U[2*1]
        x_list[2] += self.U[2*3]
        x_list[3] += self.U[2*2]

        y_list[0] += self.U[2*0+1]
        y_list[1] += self.U[2*1+1]
        y_list[2] += self.U[2*3+1]
        y_list[3] += self.U[2*2+1]
        '''
        #plt.fill(x_list+self.U[::2], y_list+self.U[1::2], c="b", alpha=0.5)
        plt.scatter(x_list+self.U[::2], y_list+self.U[1::2], c="b", alpha=0.5)
        #plt.axes().set_aspect('equal', 'datalim')

        plt.show()

    '''
    def _Bmat(self, dNdx, dNde):
        B = np.array([[dNdx[0],       0, dNdx[1],       0, dNdx[2],       0, dNdx[3],       0],
                      [0, dNde[0],       0, dNde[1],
                          0, dNde[2],       0, dNde[3]],
                      [dNde[0], dNdx[0], dNde[1], dNdx[1], dNde[2], dNdx[2], dNde[3], dNdx[3]]])
        return B
    '''

    '''
    # shape function matrix
    def _Nmat(self, xi, eta):
        N1 = 0.25*(1-xi)*(1-eta)
        N2 = 0.25*(1+xi)*(1-eta)
        N3 = 0.25*(1+xi)*(1+eta)
        N4 = 0.25*(1-xi)*(1+eta)
        N = np.array([[N1, 0, N2, 0, N3, 0, N4, 0],
                      [0, N1, 0, N2, 0, N3, 0, N4]])
        return N
    '''

    def tecplot(self):
        plt_text_header = '''
        TITLE="Topology Optimization"
        VARIABLES = "X", "Y", "UX", "UY", "FX", "FY", "rho"
        ZONE N ={0},E ={1}, DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL, VARLOCATION=([7]=CELLCENTERED)
        '''.format(self.all_node_count, self.all_element_count)
        x_list = [row[0] for row in self.node_coordinate_values]
        y_list = [row[1] for row in self.node_coordinate_values]


        with open("mesh_data/mesh.plt", "w") as plt_file:
            plt_file.write(plt_text_header)
            content = [x_list, y_list,
                       self.U[::2], self.U[1::2],
                       self.F[::2], self.F[1::2],
                       ]
            #print(content)
            np.savetxt(plt_file, content, delimiter='\n')
            np.savetxt(plt_file, self.rho.flatten(), delimiter='\n')
            np.savetxt(plt_file, self.node_connection,
                       fmt="%.0f", delimiter=',')


def main():
    start_time = time.time()
    nx = 1
    ny = 1

    vol = 1
    pnl = 3
    mesh_size = 1

    # x=0 edge fix
    fix_nodes = np.array([2*(node*(nx+1)) for node in range(ny+1)])

    # fix_nodes = np.array([0, 1, 3, 4])
    # fix_nodes = np.array([0, 1, 2, 3])

    # x=nx(*mesh_size) edge force
    F = np.zeros(2*(nx+1)*(ny+1))
    for node in range(ny+1):
        F[2*(node*(nx+1)+nx)] = 1
    # print(F)

    rho = vol*np.ones([nx, ny])
    # rho = np.random.rand(nx,ny)

    fem_obj = FEM(pnl, nx, ny, mesh_size, fix_nodes,F)
    U, l = fem_obj.fem(rho)

    print("U:", U)
    print("l:", l)
    print("elapse time [sec]:", time.time()-start_time)

    #fem_obj.plot_mesh()
    fem_obj.tecplot()


if __name__ == "__main__":
    main()
