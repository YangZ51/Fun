import pywt
import matplotlib.pyplot as plt


def plot_wavedec(input_data, wavelet_type, level, title):
    """
    Plot decomposed signals based on DWT method
    plot_dwt_decon(source[:, 'sym4', "DWT: Vibroseis Signal", 6)

    :param input_data:               Input data
    :param wavelet_type:            Wavelet type --> Ex. 'db4', 'sym4'
    :param level:                   Decomposition level
    :param title:                   Title your plot --> Ex. 'Title'
    :return:                        Decomposition plots
    """
    # DWT Decomposition
    w = pywt.Wavelet(wavelet_type)
    d = input_data
    ca = []
    cd = []
    for i in range(level):
        a, d = pywt.dwt(d, w, mode='smooth')
        ca.append(a)
        cd.append(d)

    rec_a = []
    rec_d = []

    for i, coeff in enumerate(ca):
        coeff_list = [coeff, None] + [None] * i
        rec_a.append(pywt.waverec(coeff_list, w))

    for i, coeff in enumerate(cd):
        coeff_list = [None, coeff] + [None] * i
        rec_d.append(pywt.waverec(coeff_list, w))

    fig = plt.figure(figsize=(16, 9))
    ax_main = fig.add_subplot(len(rec_a) + 1, 1, 1)
    ax_main.set_title(title)
    ax_main.plot(input_data)
    ax_main.set_xlim(0, len(input_data) - 1)

    # Plot Approximation Coefficients
    for i, y in enumerate(rec_a):
        ax = fig.add_subplot(len(rec_a) + 1, 2, 3 + i * 2)
        ax.plot(y, 'r')
        ax.set_xlim(0, len(y) - 1)
        ax.set_ylabel("A%d" % (i + 1))

    # Plot Detail Coefficients
    for i, y in enumerate(rec_d):
        ax = fig.add_subplot(len(rec_d) + 1, 2, 4 + i * 2)
        ax.plot(y, 'g')
        ax.set_xlim(0, len(y) - 1)
        ax.set_ylabel("D%d" % (i + 1))
    plt.tight_layout()
    plt.show()
