def jacobi_solver(A, b, TOL, N):
    n = len(b)
    XO = [0] * n
    #step 2
    for k in range(1, N + 1):
        x = [(1 / A[i][i]) * (b[i] - sum(A[i][j] * XO[j] for j in range(n) if j != i)) for i in range(n)]
       
        if all(abs(x[i] - XO[i]) < TOL for i in range(n)):
            print("The procedure was successful.")
            return x
        
        XO = x.copy()
    
    print("Maximum number of iterations exceeded.")
    return None

#ex
A = [[3, -2, 0],
     [-1, 5, -1],
     [0, -1, 3]]
b = [10, 8, 5]
TOL = 1e-6
N = 100

solution = jacobi_solver(A, b, TOL, N)
print("Solution:", solution)

def iterative_solver(A, b, TOL, N):
    n = len(b)
    k = 1
    XO = [0] * n
    
    while k <= N:
        #step 3
        x = [0] * n
        for i in range(n):
            sigma = sum(A[i][j] * XO[j] for j in range(i)) + sum(A[i][j] * XO[j] for j in range(i + 1, n))
            x[i] = (1 / A[i][i]) * (b[i] - sigma)
        
        #step 4
        if all(abs(x[i] - XO[i]) < TOL for i in range(n)):
            return x
        
        #step 5
        k += 1
        
        #step 6
        XO = x.copy()
    
    #step 7
    return None

#ex
A = [[3, -2, 0],
     [-1, 5, -1],
     [0, -1, 3]]
b = [10, 8, 5]
TOL = 1e-6
N = 100

solution = iterative_solver(A, b, TOL, N)
print("Solution:", solution)

def SOR_solver(A, b, TOL, N, omega):
    n = len(b)
    k = 1
    XO = [0] * n
    
    while k <= N:
        #step 3
        x = [0] * n
        for i in range(n):
            sigma = sum(A[i][j] * x[j] for j in range(i)) + sum(A[i][j] * XO[j] for j in range(i + 1, n))
            x[i] = (1 - omega) * XO[i] + (omega / A[i][i]) * (b[i] - sigma)
        
        #step 4
        if all(abs(x[i] - XO[i]) < TOL for i in range(n)):
            return x
        
        #step 5
        k += 1
        
        #step 6
        XO = x.copy()
    
    #step 7
    return None

#ex
A = [[3, -2, 0],
     [-1, 5, -1],
     [0, -1, 3]]
b = [10, 8, 5]
TOL = 1e-6
N = 100
omega = 1.2

solution = SOR_solver(A, b, TOL, N, omega)
print("Solution:", solution)

def gaussian_elimination(A, b):
    n = len(b)
    
    #forwardelim
    for i in range(n):
        for j in range(i+1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            b[j] -= factor * b[i]
    
    #backwardsub
    x = [0] * n
    for i in range(n-1, -1, -1):
        x[i] = b[i] / A[i][i]
        for j in range(i+1, n):
            x[i] -= A[i][j] / A[i][i] * x[j]
    
    return x

def iterative_refinement(A, b, TOL, N):
    n = len(b)
    k = 1
    x = gaussian_elimination(A, b)
    xx_inf = max(abs(xi) for xi in x)
    
    while k <= N:
        #step 3
        r = [b[i] - sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
        
        #step 4
        y = gaussian_elimination(A, r)
        
        #step 5
        xi = [x[i] + y[i] for i in range(n)]
        
        #step 6
        if k == 1:
            COND = max(abs(yi) for yi in y) / xx_inf
        
        #step 7
        if max(abs(xi[i] - x[i]) for i in range(n)) < TOL:
            return xi, COND
        
        #step 8
        k += 1
        
        #step 9
        x = xi
    
    #step 10
    return None, COND

#ex
A = [[3, -2, 0],
     [-1, 5, -1],
     [0, -1, 3]]
b = [10, 8, 5]
TOL = 1e-6
N = 100

solution, condition_number = iterative_refinement(A, b, TOL, N)
print("Solution:", solution)
print("Condition number:", condition_number)

