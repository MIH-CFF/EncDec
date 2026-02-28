import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from PIL import Image, ImageTk
from tkinter import messagebox as msg
import math
import numpy as np 
from collections import Counter
import os
import sys
import threading
def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS2
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def char_to_binary(c):
    return format(ord(c),'08b')

def binary_to_char(n):
    return int(n, 2)

def new_key_transformation(seed,key, target_length):

    large_key_list = list(key[::-1]+'10'+key[3:])
    seed=(seed)%256
    while len(large_key_list) < target_length:
        temp = list(key)
        for i in range(len(temp) - 1):
            if temp[i] != temp[i + 1]:
                temp[i], temp[i + 1] = temp[i + 1], temp[i]
                temp_str = ''.join(temp)
                tmp = len(temp_str)
                for j in range(tmp):
                    seed=(seed**2+i)%256
                    for k in range (0,tmp):
                        large_key_list.append('1' if temp_str[k] == temp_str[tmp-k-1] else '0' )
                    large_key_list.append('1' if seed&1 else '0' )
                    temp_str = temp_str[1:] + temp_str[0]
                    if len(large_key_list) >= target_length or i==j:
                        break
            if len(large_key_list) >= target_length:
                break
        seed=(seed**2+i)%256
        # for k in range(0,seed):
         
        temp1=char_to_binary(chr(seed))
        for i in range(8):
            large_key_list.append(temp1[i])
        temp = large_key_list[1:] + large_key_list[:1]
        key = ['0' if temp[i] == large_key_list[i] else '1' for i in range(len(temp))]
        
        if len(large_key_list) >= target_length:
            break

    return ''.join(large_key_list[:target_length])


def DNADecode(DNA, rule):
    Bits = ''
    j = 0
    if rule == 1:
        for i in range(4):
            if DNA[i] == 'A':
                Bits += '00'
            elif DNA[i] == 'G':
                Bits += '01'
            elif DNA[i] == 'C':
                Bits += '10'
            elif DNA[i] == 'T':
                Bits += '11'
            j += 2
    elif rule == 2:
        for i in range(4):
            if DNA[i] == 'A':
                Bits += '00'
            elif DNA[i] == 'C':
                Bits += '01'
            elif DNA[i] == 'G':
                Bits += '10'
            elif DNA[i] == 'T':
                Bits += '11'
            j += 2
    elif rule == 3:
        for i in range(4):
            if DNA[i] == 'T':
                Bits += '00'
            elif DNA[i] == 'G':
                Bits += '01'
            elif DNA[i] == 'C':
                Bits += '10'
            elif DNA[i] == 'A':
                Bits += '11'
            j += 2
    elif rule == 4:
        for i in range(4):
            if DNA[i] == 'T':
                Bits += '00'
            elif DNA[i] == 'C':
                Bits += '01'
            elif DNA[i] == 'G':
                Bits += '10'
            elif DNA[i] == 'A':
                Bits += '11'
            j += 2
    elif rule == 5:
        for i in range(4):
            if DNA[i] == 'C':
                Bits += '00'
            elif DNA[i] == 'T':
                Bits += '01'
            elif DNA[i] == 'A':
                Bits += '10'
            elif DNA[i] == 'G':
                Bits += '11'
            j += 2
    elif rule == 6:
        for i in range(4):
            if DNA[i] == 'C':
                Bits += '00'
            elif DNA[i] == 'A':
                Bits += '01'
            elif DNA[i] == 'T':
                Bits += '10'
            elif DNA[i] == 'G':
                Bits += '11'
            j += 2
    elif rule == 7:
        for i in range(4):
            if DNA[i] == 'G':
                Bits += '00'
            elif DNA[i] == 'T':
                Bits += '01'
            elif DNA[i] == 'A':
                Bits += '10'
            elif DNA[i] == 'C':
                Bits += '11'
            j += 2
    elif rule == 8:
        for i in range(4):
            if DNA[i] == 'G':
                Bits += '00'
            elif DNA[i] == 'A':
                Bits += '01'
            elif DNA[i] == 'T':
                Bits += '10'
            elif DNA[i] == 'C':
                Bits += '11'
            j += 2
    return Bits

def DNAEncode(bits, rule):
    DNA = ''
    j = 0
    if rule == 1:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'A'
            elif bits[i:i+2] == '01':
                DNA += 'G'
            elif bits[i:i+2] == '10':
                DNA += 'C'
            elif bits[i:i+2] == '11':
                DNA += 'T'
            j += 1
    elif rule == 2:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'A'
            elif bits[i:i+2] == '01':
                DNA += 'C'
            elif bits[i:i+2] == '10':
                DNA += 'G'
            elif bits[i:i+2] == '11':
                DNA += 'T'
            j += 1
    elif rule == 3:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'T'
            elif bits[i:i+2] == '01':
                DNA += 'G'
            elif bits[i:i+2] == '10':
                DNA += 'C'
            elif bits[i:i+2] == '11':
                DNA += 'A'
            j += 1
    elif rule == 4:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'T'
            elif bits[i:i+2] == '01':
                DNA += 'C'
            elif bits[i:i+2] == '10':
                DNA += 'G'
            elif bits[i:i+2] == '11':
                DNA += 'A'
            j += 1
    elif rule == 5:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'C'
            elif bits[i:i+2] == '01':
                DNA += 'T'
            elif bits[i:i+2] == '10':
                DNA += 'A'
            elif bits[i:i+2] == '11':
                DNA += 'G'
            j += 1
    elif rule == 6:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'C'
            elif bits[i:i+2] == '01':
                DNA += 'A'
            elif bits[i:i+2] == '10':
                DNA += 'T'
            elif bits[i:i+2] == '11':
                DNA += 'G'
            j += 1
    elif rule == 7:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'G'
            elif bits[i:i+2] == '01':
                DNA += 'T'
            elif bits[i:i+2] == '10':
                DNA += 'A'
            elif bits[i:i+2] == '11':
                DNA += 'C'
            j += 1
    elif rule == 8:
        for i in range(0, 8, 2):
            if bits[i:i+2] == '00':
                DNA += 'G'
            elif bits[i:i+2] == '01':
                DNA += 'A'
            elif bits[i:i+2] == '10':
                DNA += 'T'
            elif bits[i:i+2] == '11':
                DNA += 'C'
            j += 1
    return DNA

def DNAXOR(DNA1, DNA2):
    DNA = ''
    for i in range(4):
        if ((DNA1[i] == 'A' and DNA2[i] == 'A') or
            (DNA1[i] == 'T' and DNA2[i] == 'T') or
            (DNA1[i] == 'C' and DNA2[i] == 'C') or
            (DNA1[i] == 'G' and DNA2[i] == 'G')):
            DNA += 'A'
        elif ((DNA1[i] == 'A' and DNA2[i] == 'T') or
              (DNA1[i] == 'T' and DNA2[i] == 'A')):
            DNA += 'T'
        elif ((DNA1[i] == 'A' and DNA2[i] == 'C') or
              (DNA1[i] == 'C' and DNA2[i] == 'A')):
            DNA += 'C'
        elif ((DNA1[i] == 'A' and DNA2[i] == 'G') or
              (DNA1[i] == 'G' and DNA2[i] == 'A')):
            DNA += 'G'
        elif ((DNA1[i] == 'T' and DNA2[i] == 'C') or
              (DNA1[i] == 'C' and DNA2[i] == 'T')):
            DNA += 'G'
        elif ((DNA1[i] == 'T' and DNA2[i] == 'G') or
              (DNA1[i] == 'G' and DNA2[i] == 'T')):
            DNA += 'C'
        elif ((DNA1[i] == 'C' and DNA2[i] == 'G') or
              (DNA1[i] == 'G' and DNA2[i] == 'C')):
            DNA += 'T'
    return DNA


def encryption(ScrambleKeyBits, image_path):
    """Applies binary segments to the image's RGB channels and prints pixel values."""
    image = Image.open(image_path)
    image_array = np.array(image)
    
    #processed_image = np.copy(image)
    height, width, channels = image_array.shape
    processed_image = np.zeros(shape=(height,width,channels), dtype=np.uint8) 

    index = 0
    for i in range(height):
        for j in range(width):
            for k in range(channels):
                if index < len(ScrambleKeyBits):                    
                    rule=(index%8)+1
                    key_segment = ScrambleKeyBits[index*8:(index+1)*8]
                    original_value = image_array[i, j, k]
                    original_bin = format(original_value, '08b')

                    key_dna=DNAEncode(key_segment,rule)
                    ori_dna=DNAEncode(original_bin,rule)
                    xored_dna=DNAXOR(key_dna,ori_dna)
                    decoded_bin=DNADecode(xored_dna,rule)
                    cipher_value = int(decoded_bin, 2)

                    processed_image[i, j, k] = cipher_value
                    index += 1
                    #print(index)

    #processed_image = np.clip(processed_image, 0, 255).astype(np.uint8)

    return processed_image

def run():

    app = tk.Tk()
    global state
    global bg,fg
    bg='black'
    fg='#00FF00'
    state=0
    app.title("Image Browser")
    logo=ImageTk.PhotoImage(Image.open(resource_path("images\\logo.png")))
    app.iconphoto(True,logo)
    app.geometry("800x600+300+100")
    app.update()
    app.configure(bg=bg)
    global tk_image, tk_image2,real_img
    tk_image=None
    tk_image2=None
    real_img=None
    #functions
    def browse_image():
        global file_path
        file_path = filedialog.askopenfilename(
            filetypes=[("Image Files", "*.png *.jpg *.jpeg *.gif")])
        if file_path:
            show_image(file_path)
            
    def generate():
        if tk_image is None:
            msg.showerror("Error","Please Select image first!")
            return        
        pw = password.get()
        if len(pw) == 0:
            msg.showerror("Error","Please Input Password first!")
            return
        show_loading()
        threading.Thread(target=process_image, args=(pw,), daemon=True).start() 
    def process_image(pw):
        global tk_image2
        img = convert(pw)
        tk_image2 = ImageTk.PhotoImage(res(img,0.4,0.4))
        app.after(0, update_gui) 
        
    def update_gui():
        close_loading()
        image_label2.config(image=tk_image2, bd=2)

        con = option_var.get()

        if con == 'Encryption':
            image_1des.config(text='Plain Image')
            image_2des.config(text='Encrypted Image')
        else:
            image_2des.config(text='Decrypted Image')
            image_1des.config(text='Encrypted Image')
            
    def show_loading():
        global loading_window, progress_bar

        loading_window = tk.Toplevel(app)
        loading_window.title("Processing")
        loading_window.geometry("300x120+550+340")
        loading_window.configure(bg=bg)
        loading_window.resizable(False, False)

        label = tk.Label(
            loading_window,
            text="Processing Image...\nPlease wait",
            font=('arial',12,'bold'),
            bg=bg,
            fg=fg
        )
        label.pack(pady=10)

        progress_bar = ttk.Progressbar(
            loading_window,
            mode='indeterminate',
            length=200
        )
        progress_bar.pack(pady=10)

        progress_bar.start(10)

        loading_window.transient(app)
        loading_window.grab_set()
    
    def close_loading():
        progress_bar.stop()
        loading_window.destroy()
    
    def res(img,x,y):
        label_width = int(x * app.winfo_width())
        label_height= int(y * app.winfo_width())
        try:
            img = img.resize((label_width, label_height), Image.ANTIALIAS)
        except AttributeError:
            img = img.resize((label_width, label_height), Image.LANCZOS)
        return img  
    def download_img():
        global state
        if not state:
            msg.showerror("Error","Please select and convert image first!") 
        else:
            global real_img
            try:
                if real_img!=None:
                    download_path = filedialog.asksaveasfilename(
                    defaultextension=".png",
                    filetypes=[("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("All files", "*.*")]
        )
                    if download_path:
                        real_img.save(download_path)
                    else:
                        msg.showerror("Error","Error try again")
            except:
                msg.showerror("Error","Generate Image first")
    # Browse Button
    browse_btn = tk.Button(app, text="Browse Image",relief='groove',bd=3, command=browse_image,font=('arial','15','bold'),bg=bg,fg=fg,activebackground=bg,activeforeground=fg)
    browse_btn.place(relx=0.1, rely=0.05, relwidth=0.25, relheight=0.08)

    # Dropdown Menu (Encryption / Decryption)
    label = tk.Label(app, text="Choose operation",font=('arial','15','bold'),bg=bg,fg=fg,activebackground=bg,activeforeground=fg)
    label.place(relx=0.45, rely=0.05, relwidth=0.27,relheight=0.08)
    
    option_var = tk.StringVar(value="Encryption")
    dropdown = tk.OptionMenu(app, option_var, "Encryption", "Decryption")
    dropdown.place(relx=0.74, rely=0.05, relwidth=0.2, relheight=0.08)
    dropdown.config(relief='groove',bd=3,font=('arial','15','bold'),bg=bg,fg=fg,activebackground=bg,activeforeground=fg)
    dem1=Image.open(resource_path('images\\demo.jpg'))
    demo=ImageTk.PhotoImage(res(dem1,0.4,0.4))
    dnld1=Image.open(resource_path('images\\dw.png'))
    dnld=ImageTk.PhotoImage(res(dnld1,0.035,0.037))
    # Image Display Label
    image_label = tk.Label(app,image=demo,bg=bg)
    image_label.place(relx=0.05, rely=0.15, relwidth=0.4)
    
    image_1des = tk.Label(app,text='Original Image',font=('arial','15','bold'),bg=bg,fg=fg)
    image_1des.place(relx=0.15, rely=0.7, relwidth=0.2)
    
    image_label2 = tk.Label(app,image=demo,bg=bg)
    image_label2.place(relx=0.55, rely=0.15, relwidth=0.4)
    
    image_2des = tk.Label(app,text='Changed Image',font=('arial','15','bold'),bg=bg,fg=fg)
    image_2des.place(relx=0.65, rely=0.7, relwidth=0.2)
    download_btn = tk.Button(app,image=dnld,relief='flat',borderwidth=0, command=download_img,bg=bg,fg=fg,activebackground=bg,activeforeground=fg)
    download_btn.place(relx=0.88, rely=0.7, relwidth=0.04,relheight=0.06)
    
    entry_label=tk.Label(app,text='Input password :',font=('arial','15','bold'),bg=bg,fg=fg)
    entry_label.place(relx=0.3, rely=0.8, relwidth=0.2)
    
    password=tk.StringVar()
    entry_pass=tk.Entry(app,font=('arial','15','bold'),bg=bg,fg=fg,relief='groove',bd=3,textvariable=password,insertbackground="white")
    entry_pass.place(relx=0.55, rely=0.8, relwidth=0.3)
    
    
    
    generate_btn = tk.Button(app, text="Generate",relief='groove',bd=3, command=generate,font=('arial','15','bold'),bg=bg,fg=fg,activebackground=bg,activeforeground=fg)
    generate_btn.place(relx=0.4, rely=0.9, relwidth=0.2)
    
    def show_image(path):
        img = res(Image.open(path),0.4,0.4)
        global tk_image
        tk_image = ImageTk.PhotoImage(img)
        image_label.config(image=tk_image)
    

    def convert(pw):
        global file_path,state,real_img
        state=1
        img=Image.open(file_path)
        
        image_array = np.array(img)
        image_row, image_col, channel = image_array.shape
        target_length = image_row * image_col * channel * 8
        seed=0
        for i in pw:
            seed+=ord(i)
        binary_string = ''.join(char_to_binary(c) for c in pw)
        transformed_key= new_key_transformation(seed,binary_string, target_length)
        cipherImage =encryption(transformed_key, file_path)
        img=Image.fromarray(cipherImage)
        real_img=img
        return img
      # To keep reference
    try:
        app.mainloop()
    except KeyboardInterrupt:
        print("Application Interrupted")
    
    

if __name__ == "__main__":
    run()
    
