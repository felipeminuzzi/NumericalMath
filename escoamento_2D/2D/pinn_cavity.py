import pickle
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

tf.random.set_seed(42)
np.random.seed(42)

class PINN(tf.keras.Model):
    def __init__(self):
        super(PINN, self).__init__()
        self.hidden_layers = []
        num_hidden_layers = 8
        num_neurons_per_layer = 20

        for _ in range(num_hidden_layers):
            layer = tf.keras.layers.Dense(
                num_neurons_per_layer,
                activation='tanh',
                kernel_initializer='glorot_normal'
            )
            self.hidden_layers.append(layer)

        self.u_output = tf.keras.layers.Dense(1, activation=None)
        self.v_output = tf.keras.layers.Dense(1, activation=None)
        self.p_output = tf.keras.layers.Dense(1, activation=None)

    def call(self, inputs):
        x, y = inputs[:, 0:1], inputs[:, 1:2]
        X = tf.concat([x, y], axis=1)
        for layer in self.hidden_layers:
            X = layer(X)
        u = self.u_output(X)
        v = self.v_output(X)
        p = self.p_output(X)
        return u, v, p

def loss_fn(model, X_f, X_b, u_b, v_b):
    # Compute the PDE residuals at collocation points
    with tf.GradientTape(persistent=True) as tape:
        tape.watch(X_f)
        u, v, p = model(X_f)
        u = tf.squeeze(u)
        v = tf.squeeze(v)
        p = tf.squeeze(p)

        u_x = tape.gradient(u, X_f)[:, 0]
        u_y = tape.gradient(u, X_f)[:, 1]
        v_x = tape.gradient(v, X_f)[:, 0]
        v_y = tape.gradient(v, X_f)[:, 1]

        # Continuity equation residual
        continuity = u_x + v_y

        # Momentum equations residuals
        with tf.GradientTape(persistent=True) as tape2:
            tape2.watch(X_f)
            # Recalculate u, v within the tape2 context
            u, v, p = model(X_f)
            u = tf.squeeze(u)
            v = tf.squeeze(v)
            p = tf.squeeze(p)

            u_x = tape.gradient(u, X_f)[:, 0]
            u_y = tape.gradient(u, X_f)[:, 1]
            v_x = tape.gradient(v, X_f)[:, 0]
            v_y = tape.gradient(v, X_f)[:, 1]

            u_xx = tape2.gradient(u_x, X_f)[:, 0]
            u_yy = tape2.gradient(u_y, X_f)[:, 1]
            v_xx = tape2.gradient(v_x, X_f)[:, 0]
            v_yy = tape2.gradient(v_y, X_f)[:, 1]

        p_x = tape.gradient(p, X_f)[:, 0]
        p_y = tape.gradient(p, X_f)[:, 1]

        momentum_u = u * u_x + v * u_y + p_x - nu * (u_xx + u_yy)
        momentum_v = u * v_x + v * v_y + p_y - nu * (v_xx + v_yy)

        # PDE loss
        f_pde = tf.reduce_mean(tf.square(continuity)) + \
                tf.reduce_mean(tf.square(momentum_u)) + \
                tf.reduce_mean(tf.square(momentum_v))

    # Boundary condition loss
    u_pred_b, v_pred_b, _ = model(X_b)
    bc_loss = tf.reduce_mean(tf.square(u_pred_b - u_b)) + \
              tf.reduce_mean(tf.square(v_pred_b - v_b))

    total_loss = f_pde + bc_loss

    return total_loss

def plot_heatmap(u_k, x, y):
  plt.clf()
  plt.xlabel('x')
  plt.ylabel('y')
  plt.pcolormesh(x,y,u_k, cmap=plt.cm.jet)
  plt.colorbar()

x_min, x_max = 0.0, 1.0
y_min, y_max = 0.0, 1.0 
U_lid = 1.0  # Lid velocity  
nu = 0.01

N_f = 10000 # Number of collocation points
x_f = tf.random.uniform((N_f, 1), x_min, x_max, dtype=tf.float32)
y_f = tf.random.uniform((N_f, 1), y_min, y_max, dtype=tf.float32)
X_f_tensor = tf.concat([x_f, y_f], axis=1)    

N_b = 2000 # Number of boundary points

# Left Wall (x = 0)
x_left = x_min * tf.ones((N_b // 4, 1), dtype=tf.float32)
y_left = tf.random.uniform((N_b // 4, 1), y_min, y_max, dtype=tf.float32)

# Right Wall (x = 1)
x_right = x_max * tf.ones((N_b // 4, 1), dtype=tf.float32)
y_right = tf.random.uniform((N_b // 4, 1), y_min, y_max, dtype=tf.float32)

# Bottom Wall (y = 0)
x_bottom = tf.random.uniform((N_b // 4, 1), x_min, x_max, dtype=tf.float32)
y_bottom = y_min * tf.ones((N_b // 4, 1), dtype=tf.float32)

# Top Lid (y = 1)
x_top = tf.random.uniform((N_b // 4, 1), x_min, x_max, dtype=tf.float32)
y_top = y_max * tf.ones((N_b // 4, 1), dtype=tf.float32)

# Combine all boundary points
X_b_tensor = tf.concat([
    tf.concat([x_left, y_left], axis=1),
    tf.concat([x_right, y_right], axis=1),
    tf.concat([x_bottom, y_bottom], axis=1),
    tf.concat([x_top, y_top], axis=1)
], axis=0)

u_top = U_lid * tf.ones((N_b // 4, 1), dtype=tf.float32)
v_top = tf.zeros((N_b // 4, 1), dtype=tf.float32)

u_side = tf.zeros((3 * N_b // 4, 1), dtype=tf.float32)
v_side = tf.zeros((3 * N_b // 4, 1), dtype=tf.float32)

u_b_tensor = tf.concat([u_side, u_top], axis=0)
v_b_tensor = tf.concat([v_side, v_top], axis=0)

optimizer = tf.keras.optimizers.Adam(learning_rate=1e-3)

# Initialize the model
model = PINN()

# Compile the training step for better performance on GPU
@tf.function
def train_step():
    with tf.GradientTape() as tape:
        loss_value = loss_fn(model, X_f_tensor, X_b_tensor, u_b_tensor, v_b_tensor)
    grads = tape.gradient(loss_value, model.trainable_variables)
    optimizer.apply_gradients(zip(grads, model.trainable_variables))
    return loss_value

# Training loop
epochs = 30000
for epoch in range(epochs):
    loss_value = train_step()
    
    if epoch % 100 == 0:
        print(f"Epoch {epoch}, Loss: {loss_value.numpy():.5f}")

# Create a grid for plotting
nx, ny = 50, 50
x_plot = tf.linspace(x_min, x_max, nx)
y_plot = tf.linspace(y_min, y_max, ny)
X_grid, Y_grid = tf.meshgrid(x_plot, y_plot)
X_flat = tf.reshape(X_grid, [-1])
Y_flat = tf.reshape(Y_grid, [-1])
XY = tf.stack([X_flat, Y_flat], axis=1)

# Predict u, v at grid points
u_pred, v_pred, _ = model(XY)
u_pred = tf.reshape(u_pred, (ny, nx)).numpy()
v_pred = tf.reshape(v_pred, (ny, nx)).numpy()

save_name = ['u_pinn_pred.pkl', 'v_pinn_pred.pkl']
vels      = [u_pred, v_pred]
for i, name in enumerate(save_name):
    with open(name, 'wb') as fp:
        pickle.dump(vels[i], fp) 

# Plot streamlines
plt.figure(figsize=(8, 6))
plt.streamplot(
    X_grid.numpy(), Y_grid.numpy(),
    u_pred, v_pred, density=2, linewidth=1, arrowsize=1
)
plt.title('Predicted Velocity Field (Streamlines)')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.show()


