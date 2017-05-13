txt = file('galaxies.id.dat').read().strip()

dims, data = txt.split(':')
width, height = map(int, dims.split('x'))

print('%d %d' % (width, height))

R = 2*height-1
C = 2*width-1

out = [C*['.'] for r in range(R)]

r = 0
c = -1
for ch in data:
    steps = ord(ch) - ord('a')
    if ch != 'z': steps += 1

    for i in range(steps):
        c += 1
        if c >= C:
            r += 1
            c = 0

    if ch != 'z':
        out[r][c] = 'o'

for line in out:
    print(''.join(line))
