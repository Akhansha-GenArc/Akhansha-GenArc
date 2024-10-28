import os
import shutil
import subprocess
import time
import threading
from PyQt5.QtWidgets import QDesktopWidget
from PyQt5.QtWidgets import QApplication, QWidget, QComboBox, QLabel, QPushButton, QFileDialog, QHBoxLayout, QVBoxLayout, QMessageBox
from PyQt5.QtCore import QObject, pyqtSignal, pyqtSlot

class Signal(QObject):
    started = pyqtSignal()
    completed = pyqtSignal(float)
    error = pyqtSignal(str)

class BUSCO_GUI(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('BUSCO')
        layout = QVBoxLayout()

        # Call the center method after setting the geometry
        self.center()

        # Set up the GUI elements
        self.fasta_label = QLabel("Select a FASTA file:")
        self.fasta_button = QPushButton("Select File")
        self.fasta_file = ""
        self.fasta_path_label = QLabel()

        self.busco_label = QLabel("Select a BUSCO database for your closest species (Ex: for plants(viridiplantae))")
        self.busco_combo = QComboBox()
        with open("database.txt") as f: # Read the database names from the file
            databases = f.read().splitlines()
        self.busco_combo.addItems(databases) # Add the database names to the QComboBox

        self.output_label = QLabel("Select Output directory:")
        self.output_button = QPushButton("Select Directory")
        self.output_dir = ""
        self.output_path_label = QLabel()

        self.run_button = QPushButton("Run BUSCO")

        # Set up the layout
        fasta_layout = QHBoxLayout()
        fasta_layout.addWidget(self.fasta_label)
        fasta_layout.addWidget(self.fasta_button)
        fasta_layout.addWidget(self.fasta_path_label)

        busco_layout = QHBoxLayout()
        busco_layout.addWidget(self.busco_label)
        busco_layout.addWidget(self.busco_combo)

        output_layout = QHBoxLayout()
        output_layout.addWidget(self.output_label)
        output_layout.addWidget(self.output_button)
        output_layout.addWidget(self.output_path_label)

        layout = QVBoxLayout()
        layout.addLayout(fasta_layout)
        layout.addLayout(busco_layout)
        layout.addLayout(output_layout)
        layout.addWidget(self.run_button)
        self.setLayout(layout)

        # Connect the signals to the slots
        self.fasta_button.clicked.connect(self.select_fasta_file)
        self.busco_combo.currentIndexChanged.connect(self.update_busco_database)
        self.output_button.clicked.connect(self.select_output_dir)
        self.run_button.clicked.connect(self.run_busco)

        # Set default values for the BUSCO database
        self.busco_database = self.busco_combo.currentText()
        self.signal = Signal()
        self.signal.completed.connect(self.show_success_msg)
        self.signal.error.connect(self.show_error_msg)

    def center(self):
        """Center the window on the screen."""
        screen_rect = QDesktopWidget().availableGeometry()
        center_pos = screen_rect.center()

        self.move(center_pos.x() - self.width() // 2, center_pos.y() - self.height() // 2)


    def select_fasta_file(self):
        # Open a file dialog to select the FASTA file
        file_name, _ = QFileDialog.getOpenFileName(self, "Select File", "", "FASTA Files (*.fasta)")
        self.fasta_file = file_name
        self.fasta_path_label.setText(file_name) # Update the file path label

    def update_busco_database(self):
        # Update the BUSCO database selection
        self.busco_database = self.busco_combo.currentText()

    def select_output_dir(self):
        # Open a file dialog to select the output directory
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        self.output_dir = dir_path
        self.output_path_label.setText(dir_path) # Update the directory path label

    def run_busco(self):
        # Check if all required inputs are provided
        if not self.fasta_file:
            self.show_message_box("Missing Input", "Please select a FASTA file.")
            return
        if not self.busco_database:
            self.show_message_box("Missing Input", "Please select a BUSCO database.")
            return
        if not self.output_dir:
            self.show_message_box("Missing Output Directory", "Please select an output directory.")
            return

        start_time = time.time()  # Start time

        # Create the output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

        self.show_running_msg()

        def run_process():
            # Activate the busco_env environment and run the busco command
            activate_command = f"conda run -n busco_env busco -i {self.fasta_file} -o busco_output -l {self.busco_database} -m transcriptome --cpu 16 -f"
            process = subprocess.Popen(activate_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            if process.returncode == 0:
                # Move the output directory and downloads directory to the selected output directory
                try:
                    shutil.move("busco_output", os.path.join(self.output_dir, "busco_output"))
                    shutil.move("busco_downloads", os.path.join(self.output_dir, "busco_downloads"))
                except Exception as e:
                    self.signal.error.emit(f"Failed to move BUSCO output directories. Error message: {e}")
                    return

                # Generate the plot within the BUSCO environment
                generate_plot_command = f"conda run -n busco_env generate_plot.py -wd {self.output_dir}/busco_output"
                print(f"Running command: {generate_plot_command}")

                # Run the generate_plot command
                subprocess.run(generate_plot_command, shell=True)

                elapsed_time = time.time() - start_time  # Elapsed time
                elapsed_time_str = time.strftime("%Hh %Mm %Ss", time.gmtime(elapsed_time))  # Format elapsed time

                self.signal.completed.emit(elapsed_time)
            else:
                self.signal.error.emit("BUSCO analysis failed. Please check the console output for details.")
                print(stderr.decode())

        # Start the thread for running the process
        thread = threading.Thread(target=run_process)
        thread.start()

    @pyqtSlot(float)
    def show_success_msg(self, elapsed_time):
        elapsed_time_str = time.strftime("%Hh %Mm %Ss", time.gmtime(elapsed_time))  # Format elapsed time
        self.close()
        self.show_message_box("Success", f"BUSCO finished successfully.\nElapsed Time: {elapsed_time_str}")

    @pyqtSlot(str)
    def show_error_msg(self, error_message):
        self.close()
        self.show_message_box("BUSCO Analysis Failed", error_message)

    def show_running_msg(self):
        self.show_message_box("Running", "The process is running. Please wait until the process finishes.")

    def show_message_box(self, title, message):
        msg_box = QMessageBox(self)
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setText(message)
        msg_box.setWindowTitle(title)
        if title == "Running":
            msg_box.setStandardButtons(QMessageBox.Cancel)
        else:
            msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.show()

if __name__ == "__main__":
    app = QApplication([])
    app.setStyle('Fusion')
    gui = BUSCO_GUI()
    gui.show()
    app.exec_()