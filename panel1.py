import tensorflow as tf
import math

bst0 = tf.Variable(50,name='bst0')
tst0 = tf.Variable(5, name='tst0')
tsk0 = tf.Variable(5, name='tsk0')
bsk = tf.constant(200)
T2 = bsk

tsk = tf.clip_by_value(tsk0, 0.5, 100)
P4 = tsk
tst = tf.clip_by_value(tst0, 0.5, 100)
U14 = tst
bst = tf.clip_by_value(bst0, 0.5, 100)
AA9 = bst

mref = tf.constant(93.566)
ibchord = tf.constant(3.333333)
Esk = tf.constant(72)
H10 = Esk
Est = tf.constant(88)
H12 = Est
kk = tf.constant(2.3)
kt = tf.constant(0.18)
maxstrain = tf.constant(0.0036)
ribSpace = tf.constant(1200)
piVal = tf.constant(math.pi)

epsilonk = tf.multiply(kk, tf.pow(tf.divide(tsk, bsk), 2))
epsilont = tf.multiply(kt, tf.pow(tf.divide(tst, bst), 2))

Et = tf.add(tf.multiply(Esk, tsk),tf.multiply(Est, tf.divide(tf.multiply(bst, tst),bsk)))


Nsk = tf.multiply(Et, epsilonk)
Nst = tf.multiply(Et, epsilont)

Nxmat = tf.multiply(Et, maxstrain)

ZEAZ = tf.multiply(H12, tf.multiply(AA9, tf.multiply(U14, tf.add((tf.divide(P4,2), tf.divide(AA9,2))))))
ZEA = tf.add(tf.multiply(H12, tf.multiply(AA9, U14)), tf.multiply(H10, tf.multiply(P4, T2)))
zbar = tf.divide(ZEAZ, ZEA)
D26 = zbar

EIbar =tf.add((tf.divide(tf.multiply(H10, tf.multiply(T2,tf.pow(P4,3)))), 12), tf.add(tf.multiply(H10,tf.multiply(T2,tf.multiply(P4,tf.pow(D26,2)))), tf.add(tf.divide(tf.multiply(H12, tf.multiply(U14*tf.pow(AA9, 3))),12), tf.multiply(H12, tf.multiply(AA9, tf.multiply(U14, tf.pow((tf.add((tf.divide(AA9, 2)), tf.subtract((tf.divide(P4,2))), D26)),2)))))))
NxEuler = tf.divide(tf.multiply(tf.pow(piVal,2), EIbar), tf.multiply(tf.pow(ribSpace,2),T2))

NxSafe = tf.divide(tf.multiply(1.5, mref), tf.multiply(0.1,tf.multiply(0.4,tf.pow(ibchord,2))))

skRes = tf.divide(Nsk, NxSafe)
stRes = tf.divide(Nst, NxSafe)
matRes = tf.div(Nxmat, NxSafe)
EulerRes = tf.divide(NxEuler, NxSafe)

area = tf.add(tf.multiply(tsk, bsk), tf.multiply(tst, bst))
loss0 = tf.add(tf.divide(tf.abs(tf.subtract(1, skRes)), 0.002), tf.divide(tf.abs(tf.subtract(1, stRes)), 0.002))
loss1 = tf.cond(EulerRes >= 1.1, loss0, lambda: tf.add(loss0, tf.divide(tf.abs(tf.subtract(1.1,EulerRes)), 0.002)))
loss2 = tf.cond(matRes >= 1.1, loss1, lambda: tf.add(loss1, tf.divide(tf.abs(tf.subtract(1.1,matRes)), 0.002)))

loss = tf.add(area, loss2)
