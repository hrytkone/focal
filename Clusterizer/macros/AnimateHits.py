import imageio

nlayers = 20
images = []
for i in range(nlayers):
    filename = f"hits_{i}.png"
    images.append(imageio.imread(filename))
imageio.mimsave('hits-anim.gif', images, duration=0.5)
