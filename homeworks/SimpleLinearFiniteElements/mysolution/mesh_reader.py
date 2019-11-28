from numpy import array, empty, zeros

def mesh_reader(input_file):
  with open(input_file) as file:
    lines = file.read().split('\n')
    n_vertices = int(lines[0].split()[0])
    n_elements = int(lines[n_vertices + 1].split()[0])

    x = empty([n_vertices])
    y = empty([n_vertices])
    z = zeros([n_vertices])
    for i in range(n_vertices):
      line = lines[i + 1].split()
      x[i] = line[0]
      y[i] = line[1]
      if (len(line) > 2):
        z[i] = line[2]

    triangles = empty([n_elements, 3], dtype=int)
    for i in range(n_elements):
      line = lines[i + n_vertices + 2].split()
      triangles[i] = array(line)

  return x, y, z, triangles
