import numpy as np
import math
import cmath
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import networkx as nx


class Kuramoto_model():

    def __init__(self, graph =nx.complete_graph(5) ):   
        # number of Lorenz Oscillator
        self.set_graph(graph)
        self.set_const()
        self.init_vec =  np.random.rand(self.N_osc)* 2* np.pi 

    # set the constant
    def set_const(self, omega_cnst = 5, K_cnst=1):
        #constant must be numpy array or flout.

        #self.omega is numpy array (size = N_osc)
        if(type(omega_cnst) is float):
            self.omega =  np.ones(self.N_osc) * omega_cnst
        elif(type(omega_cnst) is int):
            self.omega =  np.ones(self.N_osc) * float(omega_cnst)
        elif(len(omega_cnst) == self.N_osc):
            self.omega = omega_cnst
        else:
            self.omega = np.random.choice(omega_cnst, self.N_osc)

        # set K
        self.K = K_cnst

        
    def set_graph(self, graph):
        self.graph = graph
        self.N_osc = len(self.graph.nodes)

    def get_neighbor_index(self):
        self.list_i_j = []* self.N_osc
        for i in range(self.N_osc):
            self.list_i_j[i] =  self.graph.succ[i]


    # set initial phase
    def set_init_state(self, init_vec):
        if (len(init_vec) != self.N_osc):
            print("length of initial phase vector must be", self.N_osc , "!")
        else:
            self.init_vec = init_vec

    # 解く
    def solve(self, dt=0.1, N_t= 100):
        # 有向グラフの場合
        if (self.graph.is_directed()):
            self.solve_directed(dt, N_t)
            return self.phases,self.Times 
        # 無向グラフの場合
        else:
            self.solve_undirected(dt, N_t)
            return self.phases, self.Times

    # 有向グラフの場合
    def dphisdt_directed(self, phis):
        dphidt = self.omega
        for i in range(0, self.N_osc):
            j_list = self.graph.pred[i]
            dphidt[i] +=  self.K / len(j_list) * np.sum(np.sin( phis[j_list] - phis[i]))
        return dphidt

    # 有向グラフの場合    
    def solve_directed(self, dt=0.1, N_t= 100):
        
        self.dt = dt
        # XYZs is numpy arr for keep date of simulation (N, 3, t)
        self.phases = np.empty((self.N_osc, N_t))
        self.Times = np.linspace(0, self.dt * N_t, N_t)
        self.phases[:, 0] = self.init_vec
        
        for ii in range(N_t - 1):
            dX1 = self.dphisdt_directed(self.phases[:, ii])
            dX2 = self.dphisdt_directed(self.phases[:, ii] + 0.5*dX1*self.dt)
            dX3 = self.dphisdt_directed(self.phases[:, ii] + 0.5*dX2*self.dt)            
            dX4 = self.dphisdt_directed(self.phases[:, ii] + dX3*self.dt)

            self.phases[:, ii+ 1] =  self.phases[:, ii] + (dX1 + 2.0 * dX2 + 2.0 * dX3 + dX4)*(self.dt/6.0)
            self.phases = self.phases % (2*np.pi)

    # 無向グラフの場合
    def dphisdt_undirected(self, phis):
        dphidt = self.omega
        for i in range(0, self.N_osc):
            j_list = self.graph.adj[i]
            dphidt[i] +=  self.K / len(j_list) * np.sum(np.sin( phis[j_list] - phis[i]))
        return dphidt

    # 無向グラフの場合    
    def solve_undirected(self, dt=0.1, N_t= 100):
        
        self.dt = dt
        # XYZs is numpy arr for keep date of simulation (N, 3, t)
        self.phases = np.empty((self.N_osc, N_t))

        self.Times = np.linspace(0, self.dt * N_t, N_t)
        self.phases[:, 0] = self.init_vec

        for ii in range(N_t - 1):
            dX1 = self.dphisdt_undirected(self.phases[:, ii])
            dX2 = self.dphisdt_undirected(self.phases[:, ii] + 0.5*dX1*self.dt)
            dX3 = self.dphisdt_undirected(self.phases[:, ii] + 0.5*dX2*self.dt)            
            dX4 = self.dphisdt_undirected(self.phases[:, ii] + dX3*self.dt)
            self.phases[:, ii+ 1] =  self.phases[:, ii] + (dX1 + 2.0 * dX2 + 2.0 * dX3 + dX4)*(self.dt/6.0)
            self.phases = self.phases % (2*np.pi)



    def get_order_parameter(self, phases=None, get_abs=True, get_phase=True):

        if (phases is None):
            m = 1 / self.N_osc * np.sum(np.exp( 1j *  self.phases), axis=0) 
            #return np.abs(m), np.angle(m)
        else:
            m = 1 / self.N_osc * np.sum(np.exp( 1j *  phases), axis=0) 
            #return np.abs(m), np.angle(m)
        if (get_abs is True and get_phase is True):
            return np.abs(m), np.angle(m)
        elif (get_abs is True):
            return np.abs(m)
        elif (get_phase is True):
            return np.angle(m)
            
    def get_angular_velocity(self):
        self.ang_vel = np.diff(self.phases, axis=1)
        self.ang_vel = np.where(self.ang_vel > 1.9*np.pi, self.ang_vel - 2*np.pi , self.ang_vel)
        self.ang_vel = np.where(self.ang_vel < - 1.9*np.pi, self.ang_vel + 2*np.pi , self.ang_vel)

        return self.ang_vel / self.dt



def main():
    sys1 = Kuramoto_model()
    sys1.set_const(np.random.normal(10, 2, size=sys1.N_osc), 0.005)
    #sys1.set_const(10, 0)

    nx.draw_networkx(sys1.graph)
    plt.show()

    Phases, Ts = sys1.solve(0.01, 2000)
    abs_op, phase_op =sys1.get_order_parameter()
    reals = np.real(np.exp(1j * Phases ))
    imgs = np.imag(np.exp(1j * Phases ))
    op_real = np.real( abs_op * np.exp(1j *phase_op) )
    op_imag = np.imag( abs_op * np.exp(1j *phase_op) )

    #print(reals.shape)
    print("t shape",Ts.shape)


    """
    for i in range(0, sys1.N_osc):
        plt.plot(Ts, Phases[i,:])
    # 表示
    plt.show()
    """



    # 描画
    fig = plt.figure(figsize=(6, 6))
    
    plt.axes().set_aspect('equal', 'datalim')
    plt.xlim(-1.2, 1.2)
    plt.xlim(-1.2, 1.2)

    def update(i):
        print(i)
        if i > 0:
            plt.cla()     
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)
        plt.plot(reals[:,i], imgs[:, i], color='b', marker=".", linestyle='None') 
        plt.plot(op_real[i], op_imag[i], color='r', marker=".", linestyle='None') 

        plt.title("t=" + str(round(Ts[i], 3)))

    anim = animation.FuncAnimation(fig, update, interval=20, frames = len(Ts))
    plt.show()

     
             
    
    plt.plot(Ts, sys1.get_order_parameter()[0])
    plt.ylim(0, 1)
    plt.title("abs of order parameter")
    plt.show()
    

if __name__ == '__main__':
    main()
