import tensorflow as tf

bsk0 = tf.Variable(initial_value=tf.random_uniform([1], 0, 10),name='bsk')
tsk0 = tf.Variable(initial_value=tf.random_uniform([1], 0, 40), name='tsk')

bsk = tf.clip_by_value(bsk0, 1, 100)
tsk = tf.clip_by_value(tsk0, 1, 100)

# Conts

# mass = tf.subtract(tf.add(tf.add(tsk, bsk), tf.multiply(tsk, bsk)), tf.divide(tsk, bsk))
mass = tsk + bsk + (tsk * bsk) - (tsk/bsk)

opt = tf.train.GradientDescentOptimizer(0.0035)
train = opt.minimize(mass)

sess = tf.Session()

init = tf.global_variables_initializer()
sess.run(init)

for step in range(4000):
    sess.run(train)
    if step % 10 == 0:
        print(step, sess.run(bsk), sess.run(tsk), sess.run(mass))