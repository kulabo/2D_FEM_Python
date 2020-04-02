import time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from pathlib import Path
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import cg
import itertools
import sys

class FEM:
    def __init__(self, nx, ny, mesh_size, fix_nodes, F):
        self.nx = nx
        self.ny = ny
        self.mesh_size = mesh_size
        self.fix_nodes = fix_nodes
        self.F = F

        self.E0 = 100

        self.all_element_count = nx*ny
        self.all_node_count = (nx+1)*(ny+1)
        self.all_vector_count = 2*(nx+1)*(ny+1)

        self.free_nodes = np.setdiff1d(
            np.arange(self.all_vector_count), fix_nodes)

        self.node_coordinate_values = np.array(
            [[x, y] for x in range(nx+1) for y in range(ny+1)])*mesh_size

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
        # Maybe two weight (for double integration) are needed
        #  but in this case, it is't necessary because they are 1.

        nu = 1/3
        # nu: poisson ratio
        # plane stress condition
        self._Dmat = self.E0 * np.array([[1, nu, 0],
                               [nu, 1, 0],
                               [0, 0, (1-nu)/2]]) / (1-nu ** 2)

    # fenite element method
    def fem(self):
        K_sp=self._Kmat_sp()

        K_free_sp = K_sp[self.free_nodes].T[
            self.free_nodes].T
        
        U = np.zeros(self.all_vector_count)
        solve_time=time.time()

        U[self.free_nodes] = spsolve(K_free_sp, sp.lil_matrix(
            self.F[self.free_nodes]).tobsr().T)
        print("solve time [sec]:", time.time()-solve_time)

        self.U = U
        return U

    # stiffness matrix
    def _Kmat(self):
        K = np.zeros([self.all_vector_count, self.all_vector_count])
        for y in tqdm(range(self.ny)):
            for x in range(self.nx):
                Ke = self._Kemat(x, y)
                top1 = (self.ny+1)*x+y
                top2 = (self.ny+1)*(x+1)+y
                elem = [2*top1, 2*top1+1, 2*top2, 2*top2+1,
                        2*top2+2, 2*top2+3, 2*top1+2, 2*top1+3]
                for index, one_elem in enumerate(elem):
                    K[elem, one_elem] += Ke[index]

        return K

    def _Kmat_sp(self):
        row = []
        col = []
        for y in tqdm(range(self.ny)):
            for x in range(self.nx):
                top1 = (self.ny+1)*x+y
                top2 = (self.ny+1)*(x+1)+y
                elem = [2*top1, 2*top1+1, 2*top2, 2*top2+1,
                        2*top2+2, 2*top2+3, 2*top1+2, 2*top1+3]
                
                # "len(elem)" is slow because it reffers elem each time.
                # Using "8" is hard-coded but fast.
                row_temp = list(itertools.chain.from_iterable([[i]*8 for i in elem]))
                col_temp = elem*len(elem)

                row.append(row_temp)
                col.append(col_temp)
        
        data = [self._Kemat(x, y) for y in tqdm(range(self.ny))
                for x in range(self.nx)]
        K_sp = sp.coo_matrix((np.array(data).flatten(), (np.array(row).flatten(), np.array(col).flatten())), shape=(
            self.all_vector_count, self.all_vector_count)).tocsr()

        return K_sp

    # stiffness matrix of one element
    def _Kemat(self, x, y):
        Ke = np.zeros((8, 8))

        for n in range(4):
            xi, eta = self.point[n]
            w = self.weight[n]
            dNdxi = self._dNdxi(eta)
            dNdeta = self._dNdeta(xi)
            J = self._Jmat(x, y, dNdxi, dNdeta)
            B = self._Bmat(J, dNdxi, dNdeta)
            Ke += w*(B.T @ self._Dmat) @ B*np.linalg.det(J)
        return Ke

    def _Jmat(self, x, y, dNdxi, dNdeta):
        '''
        node numbers for each elements
        4 --- 3
        |     |
        |     |
        1 --- 2
        '''
        #x1=self.node_coordinate_values
        x_nodes = np.array([x, x+1, x+1, x])
        y_nodes = np.array([y, y, y+1, y+1])

        # TODO use coordinate_value

        J = np.stack([dNdxi, dNdeta]) @ (
            np.stack([x_nodes, y_nodes])*self.mesh_size).T
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

        plt.scatter(x_list, y_list, c="r", alpha=0.5, label="initial")
        plt.scatter(x_list+self.U[::2], y_list +
                    self.U[1::2], c="b", alpha=0.5, label="transformed")
        plt.legend()
        plt.savefig('mesh_data/mesh.png')
        plt.show()

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
            np.savetxt(plt_file, content, delimiter='\n')
            np.savetxt(plt_file, self.rho.flatten(), delimiter='\n')
            np.savetxt(plt_file, self.node_connection,
                       fmt="%.0f", delimiter=',')


def main():
    out_mesh_data_dir=Path("mesh_data")
    if not(out_mesh_data_dir.exists()):
        out_mesh_data_dir.mkdir()
    
    start_time = time.time()
    nx = 50
    ny = 50

    mesh_size = 1

    uniformly_distributed_F = 10

    Fx = np.zeros((nx+1)*(ny+1))
    Fy = np.zeros((nx+1)*(ny+1))

    fix_x = list(range(ny+1))
    fix_y = [i for i in range(0, (nx+1)*(ny+1), ny+1)]

    Fx[-(ny+1):] = uniformly_distributed_F
    Fx[-(ny+1)] *= 0.5
    Fx[-1] *= 0.5

    fix_x = np.array(fix_x)
    fix_y = np.array(fix_y)
    fix_nodes = np.r_[2*fix_x, 2*fix_y+1]

    F = np.array([[x, y] for x, y in zip(Fx, Fy)]).flatten()

    fem_obj = FEM(nx, ny, mesh_size, fix_nodes, F)
    U = fem_obj.fem()

    print("elapse time [sec]:", time.time()-start_time)

    #fem_obj.tecplot()
    fem_obj.plot_mesh()


if __name__ == "__main__":
    main()
