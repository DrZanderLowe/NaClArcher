# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import matplotlib.pyplot as plt
import os
import sys
import tkinter as tk
import tkinter.messagebox as tkmb
import tkinter.simpledialog as tksd
import tkinter.filedialog as tkfd
import subprocess

# Main Window
class main_window(tk.Tk):
    def __init__(self, parent):
        tk.Tk.__init__(self, parent)
        self.parent=parent
        self.protocol("WM_DELETE_WINDOW", self._quit)
        self.grid()
        
        self.toolbar_frame = tk.Frame(self) 
        self.toolbar_frame.grid()
        
        self.Myfig = plt.figure(figsize=(6, 6))
        self.Canvas = FigureCanvasTkAgg(self.Myfig, master=self)
        
        self.Toolbar = NavigationToolbar2Tk(self.Canvas, self.toolbar_frame)
        self.Toolbar.update()
        self.toolbar_frame.grid(column=1,row=1,columnspan=7, sticky=tk.W+ tk.N+tk.S)
        
        self.Canvas._tkcanvas.grid(column=0,row=2, rowspan=10, columnspan=7, sticky = tk.W+tk.E+tk.N+ tk.S)
        self.result = tk.Text(self, wrap=tk.WORD, height=42, width=75)
        self.result.insert(tk.END, "")
        self.result.grid(column=0, row=2, rowspan=10, columnspan=7, sticky=tk.W + tk.E + tk.N + tk.S)
        
        self.program_1 = tk.Button(self,text = "Basic Calculation", command = self.program_1)
        self.program_2=tk.Button(self,text = "Range Calculation", command = self.program_2)
        
        self.program_1.grid(column=0, row=0,sticky=tk.E+tk.S+tk.W+tk.N)
        self.program_2.grid(column=0, row=1,sticky=tk.E+tk.S+tk.W+tk.N)

        self.QUIT = tk.Button(text = "QUIT", fg = "red", command = self._quit)
        self.QUIT.grid(column=0 , row=0, rowspan=2)
        self.QUIT.grid_remove()

        self.steps = tk.Text(self, wrap=tk.WORD,height=1, width=5)
        self.steps.insert(tk.END, "steps")
        self.steps.grid(column=1, row=0)
        
        self.temp = tk.Text(self, wrap=tk.WORD,height=1, width=10)
        self.temp.insert(tk.END, "T[K]")
        self.temp.grid(column=1, row=0)
        
        self.press = tk.Text(self, wrap=tk.WORD,height=1, width=10)
        self.press.insert(tk.END, "p[MPa]")
        self.press.grid(column=2, row=0)
        
        self.mol = tk.Text(self, wrap=tk.WORD,height=1, width=10)
        self.mol.insert(tk.END, "M[mol/kg]")
        self.mol.grid(column=3, row=0)
        
        self.calc=tk.Button(self,text = "Calculate", command = self.calculate)
        self.calc.grid(column=4, row=0,sticky=tk.E+tk.S+tk.W+tk.N)

        
        self.save=tk.Button(self,text = "Save", command = self.save_file)
        self.save.grid(column=5, row=0,sticky=tk.E+tk.S+tk.W+tk.N)
        
        self.clear_grid()
        self.result.grid()
        
        self.steps.bind("<Key>", self.clear_default_text)
        self.temp.bind("<Key>", self.clear_default_text)
        self.press.bind("<Key>", self.clear_default_text)
        self.mol.bind("<Key>", self.clear_default_text)
        
        self.steps.bind("<Return>", self.focus_next_widget)
        self.steps.bind("<Tab>", self.focus_next_widget)
        self.temp.bind("<Return>", self.focus_next_widget)
        self.temp.bind("<Tab>", self.focus_next_widget)
        self.press.bind("<Return>", self.focus_next_widget)
        self.press.bind("<Tab>", self.focus_next_widget)

    def clear_default_text(self, event):
        text = event.widget
        if (text.get("1.0", tk.END).strip() == "steps" or text.get("1.0", tk.END).strip() == "T[K]" or text.get("1.0", tk.END).strip() == "p[MPa]" or text.get("1.0", tk.END).strip() == "M[mol/kg]") : 
            text.delete("1.0", tk.END)
        

    def focus_next_widget(self, event):
        event.widget.tk_focusNext().focus()
        return "break"

    def handle_enter(self, event):
        focus_next_widget(event)
                   
    def program_1(self):
        
        self.clear_grid()

        self.result.grid()
        
        self.temp.insert(tk.END, "T[K]")
        self.temp.grid(column=2, row=0)
                
        self.press.insert(tk.END, "p[MPa]")
        self.press.grid(column=3, row=0)
        
        self.mol.insert(tk.END, "M[mol/kg]")
        self.mol.grid(column=4, row=0)
    
        self.calc.grid(column=5, row=0,sticky=tk.E+tk.S+tk.W+tk.N)
        self.save.grid(column=6, row=0,sticky=tk.E+tk.S+tk.W+tk.N)       
        


    def program_2(self):
    
        self.clear_grid()
        
        self.Canvas._tkcanvas.grid()
        self.toolbar_frame.grid()

        self.steps.insert(tk.END, "steps")
        self.steps.grid(column=1, row=0)
        
        self.temp.config(width=13)
        self.temp.insert(tk.END, "T[K]")
        self.temp.grid(column=2, row=0)
                
        self.press.insert(tk.END, "p[MPa]")
        self.press.grid(column=3, row=0)
        
        self.mol.insert(tk.END, "M[mol/kg]")
        self.mol.grid(column=4, row=0)
    
        self.calc.grid(column=5, row=0,sticky=tk.E+tk.S+tk.W+tk.N)
        self.save.grid(column=6, row=0,sticky=tk.E+tk.S+tk.W+tk.N)
    
 
    def calculate(self):
        # Provide input variables during execution

        n = self.steps.get("1.0", "end-1c")
        tempK = self.temp.get("1.0", "end-1c")
        pressure = self.press.get("1.0", "end-1c")
        molality = self.mol.get("1.0", "end-1c")
        
        n = n.replace(" ", "")
        tempK = tempK.replace(" ", "")
        pressure = pressure.replace(" ", "")
        molality = molality.replace(" ", "")
        
        try:
            if n=="":
                n=1
            n=int(n)
            if n>0:
                pass

        except ValueError:
            tkmb.showwarning("Error","Invalid input. Please input valid integer as number of steps")
            return
            
        t=tempK
        p=pressure
        m=molality
            
        tempK=[]
        pressure=[]
        molality=[]
        density_range=[]        
        
        for i in range (n):
                    
            if "-" in t:
                start, stop = t.split("-")
                start = float(start)
                stop = float(stop)
                r= (stop - start) / (n-1)
                
                tempK.append(start+i*r)
                pressure.append(p)
                molality.append(m)
                density_range.append(start+i*r)
                x_axis="Temperature [K]"
                
            if "-" in p:
                start, stop = p.split("-")
                start = float(start)
                stop = float(stop)
                r= (stop - start) / (n-1)
                
                tempK.append(t)
                pressure.append(start+i*r)
                molality.append(m)
                density_range.append(start+i*r)
                x_axis="Pressure [MPa]"
                              
            if "-" in m:
                start, stop = m.split("-")
                start = float(start)
                stop = float(stop)
                r= (stop - start) / (n-1)
                
                tempK.append(t)
                pressure.append(p)
                molality.append(start+i*r)
                density_range.append(start+i*r) 
                x_axis="Molality [mol/kg]"
                            
            if n==1:
                tempK.append(t)
                pressure.append(p)
                molality.append(m)
            

        # Check for correct input        
        try:
            for i in range (n):
                tempK[i] = float(tempK[i])
                pressure[i] = float(pressure[i])
                molality[i] = float(molality[i])

        except ValueError:
            tkmb.showwarning("Error","Invalid input. Please enter a valid number instead of T, p, M.")
            return
            
        if tempK[i]<0 or pressure[i]<0 or molality[i]<0:
            tkmb.showwarning("Error","Invalid input. Values cannot be below zero")
            return
            
        if tempK[i]==0:
            tkmb.showwarning("Error","Invalid input. Temperature cannot be equal to zero")
            return
            
        if pressure[i]==0:
            tkmb.showwarning("Error","Invalid input. Pressure cannot be equal to zero")
            return
            
        density_water=[]
        density_NaCL=[]

        
        self.data_table="Nr.\tTemperature\tPressure\tMolality\tActivity Properties\t\t\tProperties of Water\t\t\t\tProperties of NaCl Solution\t\t\t\tApparent Molar Properties\n\t[K]\t[MPa]\t[mol kg-1]\tSolute Activity Coefficient [у]\tSolvent Osmotic Activity [у]\tSolvent Activity [a]\tDensity [g cm-3]\tExpansivity [1/K]\tCompressibility [1/MPa]\tSpecific Heat Capacity [J K-1 g-1]\tDensity [g cm-3]\tExpansivity [1/K]\tCompressibility [1/MPa]\tSpecific Heat Capacity [J K-1 g-1]\tVolume [cm3 mol-1]\tExpansivity [cm3 mol-1 k-1]\tCompressibility\tRelative Enthalpy [kJ/mol]\tMolar Heat Capacity [J mol K]\n"

        os_name = os.name    
        if os_name == 'nt':
            fortran_code="fortran_code.exe"
        elif os_name == 'posix':
            fortran_code="./fortran_code"
        else:
            print(f"Unknown operating system: {os_name}")     
            return        

				# Run the compiled executable with input variables and capture output
        for i in range (n):
            result_data_values=[]
            result_process = subprocess.run(f"{fortran_code}", input=f"{tempK[i]}\n{pressure[i]}\n{molality[i]}\n", shell=True, capture_output=True, text=True)
            #print(result_process.stdout)
            # Extract numbers
            fortran_result=result_process.stdout.strip().split("\n")
            for line in fortran_result:
              if len(line)>2:
                  for j in range (len(line.split())):
                      try:
                          result_data_values.append(float(line.split() [j]))
                      except ValueError:
                          continue


            density_water.append(result_data_values[12])
            density_NaCL.append(result_data_values[11])            
                          
            self.data_table+=str(i+1)+"\t"
            for j in range(6):
              self.data_table+=str(result_data_values[j])+"\t"
            for j in range(4):
              self.data_table+=str(result_data_values[j*2+12])+"\t"
            for j in range(4):
              self.data_table+=str(result_data_values[j*2+11])+"\t"
            for j in range(5):
              self.data_table+=str(result_data_values[j+6])+"\t"
            self.data_table+="\n"
          
          
        if n==1:
            self.text_result=result_process.stdout
            self.result.grid()
        
            self.result.delete(1.0, tk.END)
            self.result.insert(tk.END, self.text_result)
            self.table_result=self.text_result
            
        else:
            self.Myfig.clf()
            ax = self.Myfig.add_subplot(111)

            ax.plot(density_range, density_water, label="Water Density")
            ax.plot(density_range, density_NaCL, label="NaCL Density")

            ax.set_xlabel(x_axis)
            ax.set_ylabel('Density [g*cm-3]')
            ax.legend()
            
            self.Canvas.draw()
                      
        # Check if the execution was successful
        if result_process.returncode == 0:
            pass
        else:
            print("Error during execution.")
            print("Error output:", result_process.stderr)
            
 
        
    def save_file(self):
        save_file=tkfd.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        with open(save_file, "w", encoding='utf-8',  newline='') as file:
            #print(self.data_table)
            file.write(self.data_table)


    def clear_grid(self):
        self.result.grid_remove()  
        self.Canvas._tkcanvas.grid_remove()  
        self.toolbar_frame.grid_remove()  

        self.steps.config(width=6)
        self.steps.delete(1.0, tk.END)
        self.steps.grid_remove()  
        
        self.temp.config(width=10)
        self.temp.delete(1.0, tk.END)
        self.temp.grid_remove()  
        
        self.press.config(width=10)        
        self.press.delete(1.0, tk.END)
        self.press.grid_remove()  
        
        self.mol.config(width=10)
        self.mol.delete(1.0, tk.END)
        self.mol.grid_remove() 
         
        self.calc.grid_remove()  
        self.save.grid_remove() 

         
    def _quit(self):
        self.destroy()
        self.quit()
        
   


# Program Start
if __name__ == "__main__":
    '''
    # Fortran source files
    fortran_files = ["nacl-2.for", "Steam.for"]
 		# Compilation command (replace gfortran with your Fortran compiler)
    compile_command = f"gfortran {' '.join(fortran_files)} -o fortran_code"
 		# Execute compilation command
    subprocess.run(compile_command, shell=True, check=True)
    '''
    app = main_window(None)
    app.title('Thermodynamic Properties Calculation')
    app.mainloop()


