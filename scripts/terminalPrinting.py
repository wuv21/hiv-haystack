from termcolor import cprint

printRed = lambda x: cprint(x, "red")
printGreen = lambda x: cprint(x, "green")
printCyan = lambda x: cprint(x, "cyan")
printBlue = lambda x: cprint(x, "blue")
printCyanOnGrey = lambda x: cprint(x, "cyan", "on_grey")

def printProgressBar(current, total, constantWord = "Progress", barLength = 30):
    percent = float(current) * 100 / total
    arrow = '-' * int(percent / 100 * barLength - 1) + '>'
    spaces = ' ' * (barLength - len(arrow))

    print('%s: [%s%s] %d %%' % (constantWord, arrow, spaces, percent), end = '\r')