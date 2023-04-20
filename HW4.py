import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, identity
from scipy.sparse.linalg import spsolve
from matplotlib.animation import FuncAnimation


## Solve the wave equation with mirroring boundary conditions
## Using centered differences in both space and time
def wavemirr(g,T,dx,dt):
    # Setup discretization
    x = np.arange(0, 1+dx, dx); t = np.arange(0, T+dt, dt)
    Nx = len(x); Nt = len(t); L = dt/dx

    # Initialize the solution matrix
    U = np.zeros((Nx, Nt))
    for k in range(1,Nx-1):
        U[k,0] = g(x[k])

    # First time step using ut(x,0)=0 and centered differences
    U[range(1,Nx-1),1] = U[range(1,Nx-1),0] + L**2/2*(U[range(2,Nx),0]-2*U[range(1,Nx-1),0]+U[range(Nx-2),0])

    # Time stepping loop using centered differences in both space and time
    for n in range(1,Nt-1):
        U[range(1,Nx-1),n+1] = 2*(1-L**2)*U[range(1,Nx-1),n] + L**2*(U[range(2,Nx),n]+U[range(Nx-2),n]) - U[range(1,Nx-1),n-1]

    return x,t,U

## Obtain the solution and create animation 'waveabs_animation.gif'
dt = 0.0001
dx = 0.001
T = 1
# Initial condition function
def g(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

x,t,U = wavemirr(g,T,dx,dt); Nx = len(x); Nt = len(t)

# Animate the solution through time
def update(n, x, U, line):
    line.set_data(x, U[:, n])
    return line,

# Set up the plot
fig, ax = plt.subplots()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('x')
ax.set_ylabel('u(x)')
ax.set_title('Absorbing Wave Equation with Centered Differences')

# Initialize the line plot
line, = ax.plot(x, U[:, 0], color='b')

# Create the animation
ani = FuncAnimation(fig, update, frames=range(0,Nt,50), fargs=(x, U, line), blit=True)

# Save the animation as a video file
ani.save('wavemirr_animation.gif', writer='pillow', fps=240)
plt.close()


## Solve the wave equation with mirroring left and reflecting right boundary
## Conditions using centered differences in both space and time
def wavemirrref(g,T,dx,dt):
    # Setup discretization
    x = np.arange(0, 1+dx, dx); t = np.arange(0, T+dt, dt)
    Nx = len(x); Nt = len(t); L = dt/dx

    # Initialize the solution matrix
    U = np.zeros((Nx, Nt))
    for k in range(1,Nx-1):
        U[k,0] = g(x[k])

    # First time step using ut(x,0)=0 and centered differences
    U[range(1,Nx-1),1] = U[range(1,Nx-1),0] + L**2/2*(U[range(2,Nx),0]-2*U[range(1,Nx-1),0]+U[range(Nx-2),0])
    # Apply boundary condition ux(1,t) = 0
    U[Nx-1,1] = U[Nx-2,1]

    # Time stepping loop using centered differences in both space and time
    for n in range(1,Nt-1):
        U[range(1,Nx-1),n+1] = 2*(1-L**2)*U[range(1,Nx-1),n] + L**2*(U[range(2,Nx),n]+U[range(Nx-2),n]) - U[range(1,Nx-1),n-1]
        # Apply boundary condition ux(1,t) = 0
        U[Nx-1,n+1] = U[Nx-2,n+1]

    return x,t,U

## Obtain the solution and create animation 'waveRref_animation.gif'
dt = 0.0001
dx = 0.001
T = 1
# Initial condition function
def g(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

x,t,U = wavemirrref(g,T,dx,dt); Nx = len(x); Nt = len(t)

# Animate the solution through time
def update(n, x, U, line):
    line.set_data(x, U[:, n])
    return line,

# Set up the plot
fig, ax = plt.subplots()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('x')
ax.set_ylabel('u(x)')
ax.set_title('Half Reflecting Wave Equation with Centered Differences')

# Initialize the line plot
line, = ax.plot(x, U[:, 0], color='b')

# Create the animation
ani = FuncAnimation(fig, update, frames=range(0,Nt,50), fargs=(x, U, line), blit=True)

# Save the animation as a video file
ani.save('wavemirrref_animation.gif', writer='pillow', fps=240)
plt.close()


# # Solve the heat equation with homogeneous Dirichlet boundary conditions using
# # Crank-Nicholson in time and centered differences in space
def heatCNCSd(g,T,dx,dt,c=1e-2,alpha=1e-6):
    # Setup discretization
    x = np.arange(0, 1+dx, dx); t = np.arange(0, T+dt, dt)
    Nx = len(x); Nt = len(t)

    # Initialize the solution matrix
    U = np.zeros((Nx, Nt))
    for k in range(1,Nx-1):
        U[k,0] = g(x[k])

    # Crank-Nicolson time-stepping
    L = c*dt/(2*dx**2)
    A = diags([-L,1+2*L,-L],[-1,0,1],shape=(Nx-2,Nx-2),format='csr')
    B = diags([L,1-2*L,L],[-1,0,1],shape=(Nx-2,Nx-2),format='csr')

    for n in range(Nt-1):
        b = B.dot(U[range(1,Nx-1),n])
        U[range(1,Nx-1),n+1] = spsolve(A,b)

        # Check for steady-state
        if np.linalg.norm(U[:, n+1] - U[:,n], np.inf) < alpha:
            tcrit = (n+1)*dt
            print(f"Steady-state reached at time {(n+1)*dt:.2f}")
            break

    return x,t[0:n],U[:,0:n],tcrit

## Plot the heat map and the total energy over the spatial domain versus time
## Display the critical time when the steady state is reached within tolerance 1e-6
dt = 0.001
dx = 0.01
T = 50
# Initial condition function
def g(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

x,t,U,tcrit = heatCNCSd(g,T,dx,dt); Nx = len(x)

A = np.zeros(len(t))
for n in range(len(t)):
    A[n] = np.trapz(U[:,n],x)

plt.imshow(np.transpose(U), extent=[0, 1, tcrit, 0], aspect='auto', cmap='hot')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("t")
plt.title("Heat equation solution using Crank-Nicolson method")
plt.show()

plt.plot(t,A,'orange')
plt.xlabel('Time')
plt.ylabel('Integrated Solution')
plt.title('Integration of Solution over Spatial Domain',y=1.05)
plt.grid(True)
plt.show()


## Solve the heat equation with homogeneous Neumann boundary conditions using
## Crank-Nicholson in time and centered differences in space
def heatCNCSn(g,T,dx,dt,c=1e-2,alpha=1e-6):
    # Setup discretization
    x = np.arange(0, 1+dx, dx); t = np.arange(0, T+dt, dt)
    Nx = len(x); Nt = len(t)

    # Initialize the solution matrix
    U = np.zeros((Nx, Nt))
    for k in range(1,Nx-1):
        U[k,0] = g(x[k])

    # Crank-Nicolson time-stepping
    L = c*dt/(2*dx**2)
    A = diags([-L,1+2*L,-L],[-1,0,1],shape=(Nx,Nx),format='csr')
    B = diags([L,1-2*L,L],[-1,0,1],shape=(Nx,Nx),format='csr')

    # Neumann boundary conditions
    A[0,1] = -2*L; B[0,1] = 2*L
    A[-1,-2] = -2*L; B[-1,-2] = 2*L

    for n in range(Nt-1):
        b = B.dot(U[:,n])
        U[:,n+1] = spsolve(A,b)

        # Check for steady-state
        if np.linalg.norm(U[:, n+1] - U[:,n], np.inf) < alpha:
            tcrit = (n+1)*dt
            print(f"Steady-state reached at time {(n+1)*dt:.2f}")
            break

    return x,t[0:n],U[:,0:n],tcrit

## Plot the heat map and the total energy over the spatial domain versus time
## Display the critical time when the steady state is reached within tolerance 1e-6
dt = 0.001
dx = 0.01
T = 50
# Initial condition function
def g(x):
    return float(np.piecewise(x,[.4<x<=.5, .5<x<.6, 0<=x<=.4 or .6<=x<=1],[lambda z: 10*z-4, lambda z:6-10*z, 0]))

x,t,U,tcrit = heatCNCSn(g,T,dx,dt); Nx = len(x)

A = np.zeros(len(t))
for n in range(len(t)):
    A[n] = np.trapz(U[:,n],x)

plt.imshow(np.transpose(U), extent=[0, 1, tcrit, 0], aspect='auto', cmap='hot')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("t")
plt.title("Heat equation solution using Crank-Nicolson method")
plt.show()

plt.plot(t,A,'orange')
plt.ylim(0,.11)
plt.xlabel('Time')
plt.ylabel('Integrated Solution')
plt.title('Integration of Solution over Spatial Domain',y=1.05)
plt.grid(True)
plt.show()