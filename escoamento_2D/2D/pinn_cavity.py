import numpy as np
import tensorflow as tf
from tensorflow.keras import layers
import matplotlib.pyplot as plt


#define PINN model
class PINN(tf.keras.Model):
    def __init__(self, layers_sizes):
        super(PINN, self).__init__()
        self.hidden = [layers.Dense(size, activation='tanh',
                                    kernel_initializer=tf.keras.initializers.GlorotNormal()) 
                       for size in layers_sizes]
        self.out = layers.Dense(3, kernel_initializer='zeros')

    def call(self, x):
        z = x
        for layer in self.hidden:
            z = layer(z)
        return self.out(z)
        
#calculate the residuals:
def residuals(model, x, y, t, nu):
    with tf.GradientTape(persistent=True) as tape:
        tape.watch([x, y, t])
        inputs = tf.stack([x, y, t], axis=1)
        outputs = model(inputs)
        u = outputs[:, 0]
        v = outputs[:, 1]
        p = outputs[:, 2]

        u_x = tape.gradient(u, x)
        u_y = tape.gradient(u, y)
        u_t = tape.gradient(u, t)

        v_x = tape.gradient(v, x)
        v_y = tape.gradient(v, y)
        v_t = tape.gradient(v, t)

        p_x = tape.gradient(p, x)
        p_y = tape.gradient(p, y)

    u_xx = tape.gradient(u_x, x)
    u_yy = tape.gradient(u_y, y)
    v_xx = tape.gradient(v_x, x)
    v_yy = tape.gradient(v_y, y)

    # Navier-Stokes residuals
    f_u = u_t + u * u_x + v * u_y + p_x - nu * (u_xx + u_yy)
    f_v = v_t + u * v_x + v * v_y + p_y - nu * (v_xx + v_yy)

    # Continuity residual
    f_cont = u_x + v_y

    return f_u, f_v, f_cont

#defining loss function
def loss_fn(model, x, y, t, nu):
    f_u, f_v, f_cont = residuals(model, x, y, t, nu)

    # PDE residual loss
    loss_pde = tf.reduce_mean(tf.square(f_u)) + tf.reduce_mean(tf.square(f_v)) + tf.reduce_mean(tf.square(f_cont))

    # Boundary conditions
    inputs = tf.stack([x, y, t], axis=1)
    outputs = model(inputs)
    u = outputs[:, 0]
    v = outputs[:, 1]

    # Lid condition: y=1, u=1, v=0
    lid = tf.where(tf.equal(y, 1.0))
    loss_bc_top = tf.reduce_mean(tf.square(tf.gather_nd(u, lid) - 1.0)) + tf.reduce_mean(tf.square(tf.gather_nd(v, lid)))

    # Other walls: u=0, v=0
    walls = tf.where(tf.not_equal(y, 1.0))
    loss_bc_walls = tf.reduce_mean(tf.square(tf.gather_nd(u, walls))) + tf.reduce_mean(tf.square(tf.gather_nd(v, walls)))

    return loss_pde + loss_bc_top + loss_bc_walls

#training
def train(model, x, y, t, nu, epochs, lr):
    optimizer = tf.keras.optimizers.Adam(learning_rate=lr)

    for epoch in range(epochs):
        with tf.GradientTape() as tape:
            loss = loss_fn(model, x, y, t, nu)
        grads = tape.gradient(loss, model.trainable_variables)
        
        grads = [tf.clip_by_norm(g, 1.0) for g in grads]
        
        optimizer.apply_gradients(zip(grads, model.trainable_variables))

        if epoch % 500 == 0:
            print(f"Epoch {epoch}: Loss = {loss.numpy()}")


#running the model:
# Define domain points
N = 5000
x = tf.random.uniform((N,), 0.0, 1.0)
y = tf.random.uniform((N,), 0.0, 1.0)
t = tf.random.uniform((N,), 0.0, 1.0)  # Time domain [0,1]

nu = 0.01  # Kinematic viscosity

# Initialize model
model = PINN([50, 50, 50])

# Train
train(model, x, y, t, nu, epochs=5000, lr=1e-4)

#visualizae:
x_plot = np.linspace(0, 1, 30)
y_plot = np.linspace(0, 1, 30)
X, Y = np.meshgrid(x_plot, y_plot)
T = np.full_like(X.flatten(), 0.5)  # mid-time slice

XYT = np.vstack([X.flatten(), Y.flatten(), T]).T
pred = model(tf.convert_to_tensor(XYT, dtype=tf.float32))

U = pred[:, 0].numpy().reshape(X.shape)
V = pred[:, 1].numpy().reshape(Y.shape)

plt.figure(figsize=(6,6))
plt.quiver(X, Y, U, V)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Velocity field at t=0.5')
plt.show()


