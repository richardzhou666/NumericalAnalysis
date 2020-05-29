class heat:
    # This class solves the classic 1d heat equation ut = uxx
    # Initial condition u(0,x) = e^(-100x^2)
    # Boundary condition u(t,-1) = u(t,1) = 0
    # Inputs
    def __init__(self, k, h):
        self.k = k # time step
        self.h = h # mesh size
        self.m = int(2 / self.h) # horizontal grid
        self.n = int(1 / self.k) # vertical grid
        self.x = np.linspace(-1, 1, self.m)
        self.alpha = self.k / self.h ** 2 

    def forwardEuler(self):
        grid = np.zeros((self.n, self.m))
        grid[0, :] = np.exp(-100*self.x**2)
        for i in range(1, self.n):
            for j in range(1, self.m-1):
                grid[i, j] = grid[i-1, j] + self.alpha*(grid[i-1, j + 1] - 2 * grid[i - 1, j] + grid[i - 1, j -1])
        return grid

    def backwardEuler(self):
        grid = np.zeros((self.n, self.m))
        grid[0, :] = np.exp(-100*self.x**2)
        k = np.array([-self.alpha*np.ones(self.m-1),(1 + 2 * self.alpha)*np.ones(self.m),-self.alpha * np.ones(self.m-1)])
        offset = [-1, 0, 1]
        A = diags(k, offset).toarray()
        for i in range(1, self.n):
            grid[i, :] = np.linalg.solve(A, grid[i-1,:])
        return grid

    def crankNicolson(self):
        grid = np.zeros((self.n, self.m))
        grid[0, :] = np.exp(-100*self.x**2)
        offset = [-1,0,1]
        k = np.array([-self.alpha / 2 *np.ones(self.m-1),(1 + self.alpha)*np.ones(self.m),-self.alpha / 2 * np.ones(self.m-1)])
        A = diags(k,offset).toarray()
        l = np.array([self.alpha / 2 *np.ones(self.m-1),(1 - self.alpha)*np.ones(self.m), self.alpha / 2 * np.ones(self.m-1)])
        B = diags(l,offset).toarray()
        for i in range(1, self.n):
            grid[i, :] = np.linalg.solve(A, np.dot(B,grid[i-1,:]))
        return grid
    
    # Setup subplots before calling plot
    # i.e. fig, ax = plt.subplots(2, 2, figsize = (10, 8))
    # Assign row to i
    def plot(self, i):
        ax[i, 0].plot(self.x, self.forwardEuler()[1, :])
        ax[i, 0].set_title("Forward Euler k = {}".format(self.k))
        ax[i, 1].plot(self.x, self.backwardEuler()[1, :])
        ax[i, 1].set_title("Backward Euler k = {}".format(self.k))
        ax[i, 2].plot(self.x, self.crankNicolson()[1, :])
        ax[i, 2].set_title("Crank Nicolson k = {}".format(self.k))