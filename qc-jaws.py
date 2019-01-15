import os
import sys
sys.path.append('C:\Program Files\GDCM 2.8\lib')
import pydicom
import subprocess
# import gdcm
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy import signal
from scipy.signal import find_peaks,peak_prominences,peak_widths

# np.set_printoptions(threshold=np.nan)

#*************FILTERING EXAMPLE**************************************
# def butter_lowpass(cutoff, fs, order=5):
#     nyq = 0.5 * fs
#     normal_cutoff = cutoff / nyq
#     b, a = butter(order, normal_cutoff, btype='low', analog=False)
#     return b, a
#
# def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
#     b, a = butter_lowpass(cutoff, fs, order=order)
#     y = filtfilt(b, a, data)
#     return y

    # cutoff = 500
    # fs = 70000
    # amp_filt = butter_lowpass_filtfilt(amplitude, cutoff, fs)

    # N=100
    # amp_filt=np.convolve(amplitude, np.ones((N,)) / N, mode='full')
#*************FILTERING EXAMPLE**************************************

def running_mean(x, N):
    out = np.zeros_like(x, dtype=np.float64)
    dim_len = x.shape[0]
    for i in range(dim_len):
        if N%2 == 0:
            a, b = i - (N-1)//2, i + (N-1)//2 + 2
        else:
            a, b = i - (N-1)//2, i + (N-1)//2 + 1

        #cap indices to min and max indices
        a = max(0, a)
        b = min(dim_len, b)
        out[i] = np.mean(x[a:b])
    return out


#axial visualization and scrolling
def multi_slice_viewer(volume,dx,dy):
    # remove_keymap_conflicts({'j', 'k'})
    fig, ax = plt.subplots()
    ax.volume = volume
    ax.index = volume.shape[2] // 2
    print(ax.index)
    extent = (0, 0 + (volume.shape[1] * dx),
              0, 0 + (volume.shape[0] * dy))
    ax.imshow(volume[:,:,ax.index],extent=extent)
    ax.set_xlabel('x distance [mm]')
    ax.set_ylabel('y distance [mm]')
    ax.set_title("slice="+ str(ax.index))
    fig.suptitle('Axial view', fontsize=16)
    fig.canvas.mpl_connect('key_press_event', process_key_axial)

def process_key_axial(event):
    fig = event.canvas.figure
    ax = fig.axes[0]
    if event.key == 'j':
        previous_slice_axial(ax)
    elif event.key == 'k':
        next_slice_axial(ax)
    ax.set_title("slice="+str(ax.index))
    fig.canvas.draw()

def previous_slice_axial(ax):
    volume = ax.volume
    ax.index = (ax.index - 1) % volume.shape[2]  # wrap around using %
    print(ax.index,volume.shape[2])
    ax.images[0].set_array(volume[:,:,ax.index])

def next_slice_axial(ax):
    volume = ax.volume
    ax.index = (ax.index + 1) % volume.shape[2]
    print(ax.index,volume.shape[2])
    ax.images[0].set_array(volume[:,:,ax.index])





def minimize_junction(amplitude, peaks,dx):
    for j in range(0, amplitude.shape[1]-1):
        cumsum_prev = 1e6
        # if j==1:
        #     exit(0)
        print('peak=',j,peaks[j])
        amp_base_res = amplitude[:,j]
        amp_overlay_res = amplitude[:,j+1]
        for i in range(0,-38,-1):
            x = np.linspace(0, 0 + (len(amp_base_res) * dx), len(amplitude), endpoint=False)  # definition of the distance axis
            amp_overlay_res_roll=np.roll(amp_overlay_res,i)

            #amplitude is the vector to analyze +-500 samples from the center
            amp_tot = ( amp_base_res[peaks[j] - 1000:peaks[j] + 1000] + amp_overlay_res_roll[peaks[j] - 1000:peaks[j] + 1000] )  #divided by 2 to normalize
            xsel = x[peaks[j] - 1000:peaks[j] + 1000]

            amp_filt = running_mean(amp_tot,281)
            cumsum=np.sum(np.abs(amp_tot - amp_filt))
            # print(i,abs(i)*dx,cumsum,cumsum_prev)

            # plt.figure(figsize=(10, 6))
            # plt.plot(amp_tot)
            # plt.plot(amp_filt)
            # plt.xlabel('x distance [mm]')
            # plt.ylabel('amplitude')
            # plt.title(str(i))


            if cumsum>cumsum_prev: # then we went too far
                print('peak=',j,'i=', i+1, "dx=", dx, "delta=", abs(i+1) * dx, "cumsum=", cumsum, "cumsum_prev=", cumsum_prev,'<-final')
                plt.figure(figsize=(10, 6))
                # plt.plot(-1*(amp_base_res+amp_overlay_res_roll))
                plt.plot(-1 * (amp_prev)) #here we multiply by -1 just to flip the graph
                plt.plot(-1*amp_filt_prev)
                plt.xlabel('x distance [mm]')
                plt.ylabel('amplitude')
                plt.title('final')
                plt.show()
                # print(amp_prev)
                break
            else:
                print('i=',i,"dx=", dx, "delta=", abs(i) * dx, "cumsum=", cumsum, "cumsum_prev=", cumsum_prev)
                amp_prev=amp_tot
                amp_filt_prev=amp_filt
                cumsum_prev=cumsum
        plt.show()












def peak_find(amp_peak):
    peaks,_ = find_peaks(amp_peak,prominence=15,height=(14000,15500),width=(50,250))
    peaks_w = peak_widths(amp_peak,peaks,rel_height=0.5)
    print('peaks_width=',peaks_w)
    print('peak prominence=',peak_prominences(amp_peak,peaks))
    plt.figure()
    plt.plot(amp_peak)
    plt.plot(peaks,amp_peak[peaks],"x")
    plt.hlines(*peaks_w[1:], color="C2")
    plt.show()
    print('peaks',peaks)
    # exit(0)
    return peaks















#this subroutine will merge the two jaws into a single image and display a graph of the overlap
def merge_view_vert(volume,dx,dy):


    #creating merged volume
    merge_vol = np.zeros((volume.shape[0], volume.shape[1]))

    #creating vector for processing along cols (one row)
    amplitude = np.zeros((volume.shape[1],volume.shape[2]))  # 1 if it is vertical 0 if the bars are horizontal

    x = np.linspace(0,0+(volume.shape[1]*dx),volume.shape[1],endpoint=False)#definition of the distance axis
    # x = np.arange(0,)#definition of the distance axis

    #merging the two images together
    ampl_resamp = np.zeros(((volume.shape[1])*10, volume.shape[2]))
    amp_peak = np.zeros((volume.shape[1]) * 10)

    for slice in tqdm(range(0,volume.shape[2])):
        merge_vol = merge_vol + volume[:, :, slice]
        amplitude[:,slice]=volume[int(volume.shape[0] / 2), :, slice]
        ampl_resamp[:,slice]=signal.resample(amplitude[:,slice],int(len(amplitude))*10) #resampling the amplitude vector
        amp_peak=amp_peak+ampl_resamp[:,slice]/volume.shape[2]





    fig, ax = plt.subplots(nrows=2,squeeze=True,sharex=True)

    extent = (0, 0 + (volume.shape[1] * dx),
              0, 0 + (volume.shape[0] * dy))

    print(len(x),merge_vol.shape)
    ax[0].imshow(merge_vol,extent=extent,aspect='auto')
    # ax[0].set_aspect('equal', 'box')
    ax[0].set_xlabel('x distance [mm]')
    ax[0].set_ylabel('y distance [mm]')

    ax[1].plot(x,amplitude)
    ax[1].set_xlabel('x distance [mm]')
    ax[1].set_ylabel('amplitude')
    # ax.set_title("slice=" + str(ax.index))
    fig.suptitle('Merged volume', fontsize=16)







    peaks = peak_find(amp_peak)
    minimize_junction(ampl_resamp,peaks,dx/10)
































#this subroutine will merge the 4 jaws and analyze the two upper and two lower pairs
#each jaw is 6cm in length (60mm)
def merge_view_horz(volume,dx,dy):


    # creating merged volume
    merge_vol = np.zeros((volume.shape[0], volume.shape[1]))

    # creating vector for processing along cols (one row)
    amplitude = np.zeros((volume.shape[0],volume.shape[2]))  # 1 if it is vertical 0 if the bars are horizontal

    x = np.linspace(0, 0 + (volume.shape[0] * dy), volume.shape[0], endpoint=False)  # definition of the distance axis
    # x = np.arange(0,) #definition of the distance axis

    # merging the two images together
    ampl_resamp = np.zeros(((volume.shape[0]) * 10, volume.shape[2]))
    amp_peak = np.zeros((volume.shape[0]) * 10)


    for slice in tqdm(range(0, volume.shape[2])):
        merge_vol = merge_vol + volume[:, :, slice]
        amplitude[:, slice] = volume[:,int(volume.shape[1] / 2), slice]
        ampl_resamp[:, slice] = signal.resample(amplitude[:, slice], int(len(amplitude)) * 10)  # resampling the amplitude vector
        amp_peak = amp_peak + ampl_resamp[:, slice] / volume.shape[2]


    fig, ax = plt.subplots(nrows=2, squeeze=True, sharey=True)


    extent = (0, 0 + (volume.shape[1] * dx),
              0, 0 + (volume.shape[0] * dy))

    print(len(x), merge_vol.shape)
    ax[0].imshow(merge_vol, extent=extent, aspect='auto')
    # ax[0].set_aspect('equal', 'box')
    ax[0].set_xlabel('x distance [mm]')
    ax[0].set_ylabel('y distance [mm]')

    ax[1].plot(amplitude,x)
    # ax[1].set_xlabel('x distance [mm]')
    # ax[1].set_ylabel('amplitude')
    ax[1].set_xlabel('amplitude')
    ax[1].set_ylabel('y distance [mm]')
    # ax.set_title("slice=" + str(ax.index))
    fig.suptitle('Merged volume', fontsize=16)





    peaks = peak_find(amp_peak)
    minimize_junction(ampl_resamp,peaks, dx / 10)











































































#the data loaded is not in HU we need to rescale it
# def get_pixels_hu(scans):
#     image = np.stack([s.pixel_array for s in scans])
#     # Convert to int16 (from sometimes int16),
#     # should be possible as values should always be low enough (<32k)
#     image = image.astype(np.int16)
#
#     # Set outside-of-scan pixels to 1
#     # The intercept is usually -1024, so air is approximately 0
#     image[image == -2000] = 0
#
#     # Convert to Hounsfield units (HU)
#     intercept = scans[0].RescaleIntercept
#     slope = scans[0].RescaleSlope
#
#     if slope != 1:
#         image = slope * image.astype(np.float64)
#         image = image.astype(np.int16)
#
#     image += np.int16(intercept)
#
#     return np.array(image, dtype=np.int16)

# interpolate the dataset so the spacing is equal in every direction
# def resample(image, scan, new_spacing=[1, 1, 1]): #this will resample the data to 1mmx1mmx1mm
#     # Determine current pixel spacing
#     spacing = map(float, ([scan[0].SliceThickness] + scan[0].PixelSpacing))
#     spacing = np.array(list(spacing))
#
#     resize_factor = spacing / new_spacing
#     new_real_shape = image.shape * resize_factor
#     new_shape = np.round(new_real_shape)
#     real_resize_factor = new_shape / image.shape
#     new_spacing = spacing / real_resize_factor
#
#     image = scipy.ndimage.interpolation.zoom(image, real_resize_factor)
#
#     return image, new_spacing


def read_dicom3D(dirname,poption):
    slice=0
    # lstFilesDCM = [] #empty list to store dicom files
    for subdir, dirs, files in os.walk(dirname):
        k = 0
        for file in tqdm(sorted(files)):
            if poption.startswith(('y', 'yeah', 'yes')):
                subprocess.call(["gdcmconv", "-w", dirname + file, os.path.splitext(dirname + file)[0] + "_decomp" + ".dcm"])
                dataset = pydicom.dcmread(os.path.splitext(dirname + file)[0] + "_decomp" + ".dcm")
                if k==0:
                    ArrayDicom = np.zeros((dataset.Rows, dataset.Columns), dtype=dataset.pixel_array.dtype)
                    ArrayDicom = np.dstack((ArrayDicom, dataset.pixel_array))
                    # print("slice thickness [mm]=",dataset.SliceThickness)
                    print("pixel spacing row [mm]=", dataset.ImagePlanePixelSpacing[0])
                    print("pixel spacing col [mm]=", dataset.ImagePlanePixelSpacing[1])
                    dx=dataset.ImagePlanePixelSpacing[0]
                    dy=dataset.ImagePlanePixelSpacing[1]
                else:
                    ArrayDicom=np.dstack((ArrayDicom, dataset.pixel_array))
            elif poption.startswith(('n', 'no', 'nope')):
                dataset = pydicom.dcmread(dirname + file)
                if k==0:
                    ArrayDicom = np.zeros((dataset.Rows,dataset.Columns,0), dtype=dataset.pixel_array.dtype)
                    ArrayDicom = np.dstack((ArrayDicom, dataset.pixel_array))
                    # print("slice thickness [mm]=", dataset.SliceThickness)
                    print("pixel spacing row [mm]=", dataset.ImagePlanePixelSpacing[0])
                    print("pixel spacing col [mm]=", dataset.ImagePlanePixelSpacing[1])
                    dx=dataset.ImagePlanePixelSpacing[0]
                    dy=dataset.ImagePlanePixelSpacing[1]
                else:
                    ArrayDicom = np.dstack((ArrayDicom, dataset.pixel_array))
            print(k)
            k=k+1



    # we have loaded the image but the spacing is different in the x,y and z directions
    #spacing = map(float, ([scan[0].SliceThickness] + scan[0].PixelSpacing))
    #spacing = np.array(list(spacing))
    print(type(ArrayDicom),np.shape(ArrayDicom))

    multi_slice_viewer(ArrayDicom,dx,dy)
    if k==2:
        merge_view_vert(ArrayDicom, dx, dy)
    else:
        merge_view_horz(ArrayDicom, dx, dy)

    plt.show(block=True)





    # Normal mode:
    print()
    print("Directory folder.........:", dirname)
    print("Storage type.....:", dataset.SOPClassUID)
    print()

    pat_name = dataset.PatientName
    display_name = pat_name.family_name + ", " + pat_name.given_name
    print("Patient's name...:", display_name)
    print("Patient id.......:", dataset.PatientID)
    print("Modality.........:", dataset.Modality)
    print("Study Date.......:", dataset.StudyDate)
    print("Gantry angle......", dataset.GantryAngle)
    #
    # if 'PixelData' in dataset:
    #     rows = int(dataset.Rows)
    #     cols = int(dataset.Columns)
    #     print("Image size.......: {rows:d} x {cols:d}, {size:d} bytes".format(
    #         rows=rows, cols=cols, size=len(dataset.PixelData)))
    #     if 'PixelSpacing' in dataset:
    #         print("Pixel spacing....:", dataset.PixelSpacing)
    #
    # # use .get() if not sure the item exists, and want a default value if missing
    # print("Slice location...:", dataset.get('SliceLocation', "(missing)"))



    # print(type(dataset))
    # print(dataset.pixel_array)
    # # exit(0)
    # # plot the image using matplotlib
    # plt.imshow(dataset.pixel_array, cmap=plt.cm.bone)
    # plt.xlabel('X Distance from reference point')
    # plt.ylabel('Y Distance from reference point')
    # plt.show(block=True)



try:
    dirname = str(sys.argv[1])
    print(dirname)
except:
    print('Please enter a valid filename')
    print("Use the following command to run this script")
    print("python test_pydicom3D.py \"[dirname]\"")


while True:  # example of infinite loops using try and except to catch only numbers
    line = input('Are the files compressed [yes(y)/no(n)]> ')
    try:
        ##        if line == 'done':
        ##            break
        poption = str(line.lower())
        if poption.startswith(('y', 'yeah', 'yes','n','no','nope')):
            break

    except:
        print('Please enter a valid option:')

read_dicom3D(dirname,poption)

