import matplotlib.pyplot as plt
import pywt


def plot_waverec(input_data, Fs, w, level):
    """
    Plot denoised signals based on DWT
    plot_waverec(input_data, 250, 'db4', 6)

    :param input_data:  Array-like data
    :param Fs:          Sampling frequency
    :param w:           Waveform
    :param title:       Plot title
    :param level:       Level of the DWT decomposition
    :return:            Denoised Signal
    """

    time_end = len(input_data)/Fs
    print("Time series length is : ", time_end, "sec")
    t = np.linspace(0, time_end, len(input_data))
    w = pywt.Wavelet(w)
    d = input_data
    d_coeffs = pywt.wavedec(d, w, mode='smooth', level=level)
    rec_coeffs = []

    for i, coeffs in enumerate(d_coeffs):
        rc = pywt.threshold_firm(coeffs, value_low=0.3, value_high=1.3)
        rec_coeffs.append(rc)

    recon_data = pywt.waverec(rec_coeffs, w, mode='smooth')

    fig, axs = plt.subplots(3, 1, figsize=(13, 9))
    axs[0].plot(t, input_data, 'k')
    axs[0].set_title('Original Signal')
    axs[1].plot(t, recon_data, 'k')
    axs[1].set_title('Denoised Signal')
    axs[2].plot(t, recon_data - input_data, 'k')
    axs[2].set_title('Difference')

    plt.tight_layout()
    plt.show()
