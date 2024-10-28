import sys
from PyQt5.QtWidgets import (
    QApplication, QDialog, QGridLayout, QLabel, QPushButton,
    QFileDialog, QMessageBox, QVBoxLayout, QTextEdit, QHBoxLayout, QSizePolicy
)
from PyQt5.QtCore import QElapsedTimer, Qt, QObject, pyqtSignal, pyqtSlot
import subprocess
import threading
import time
import shutil
from datetime import datetime
import os


class Signal(QObject):
    started = pyqtSignal()
    completed = pyqtSignal(float)
    error = pyqtSignal(str)


class PopupWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.forward_files = []
        self.reverse_files = []
        self.output_folder = ""
        self.process = None  # Store the subprocess instance
        self.signal = Signal()  # Signal object for communicating between threads
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Trinity')
        self.setFixedSize(600, 400)  # Set the window size
        layout = QVBoxLayout()

                # Forward file selection
        forward_layout = QVBoxLayout()
        forward_layout.addWidget(QLabel('Select Forward Files:'))

        self.forward_paths = QTextEdit()
        self.forward_paths.setReadOnly(True)
        forward_layout.addWidget(self.forward_paths)

        forward_button_layout = QHBoxLayout()
        forward_select_button = QPushButton('Select', self)
        forward_select_button.clicked.connect(self.select_forward_files)
        forward_button_layout.addWidget(forward_select_button)

        forward_add_button = QPushButton('Add', self)
        forward_add_button.clicked.connect(self.add_forward_files)
        forward_button_layout.addWidget(forward_add_button)

        forward_delete_button = QPushButton('Delete', self)
        forward_delete_button.clicked.connect(self.delete_forward_file)
        forward_button_layout.addWidget(forward_delete_button)

        forward_clear_button = QPushButton('Clear', self)
        forward_clear_button.clicked.connect(self.clear_forward_files)
        forward_button_layout.addWidget(forward_clear_button)

        forward_layout.addLayout(forward_button_layout)
        layout.addLayout(forward_layout)

        # Reverse file selection
        reverse_layout = QVBoxLayout()
        reverse_layout.addWidget(QLabel('Select Reverse Files:'))

        self.reverse_paths = QTextEdit()
        self.reverse_paths.setReadOnly(True)
        reverse_layout.addWidget(self.reverse_paths)

        reverse_button_layout = QHBoxLayout()
        reverse_select_button = QPushButton('Select', self)
        reverse_select_button.clicked.connect(self.select_reverse_files)
        reverse_button_layout.addWidget(reverse_select_button)

        reverse_add_button = QPushButton('Add', self)
        reverse_add_button.clicked.connect(self.add_reverse_files)
        reverse_button_layout.addWidget(reverse_add_button)

        reverse_delete_button = QPushButton('Delete', self)
        reverse_delete_button.clicked.connect(self.delete_reverse_file)
        reverse_button_layout.addWidget(reverse_delete_button)

        reverse_clear_button = QPushButton('Clear', self)
        reverse_clear_button.clicked.connect(self.clear_reverse_files)
        reverse_button_layout.addWidget(reverse_clear_button)

        reverse_layout.addLayout(reverse_button_layout)
        layout.addLayout(reverse_layout)


        # Output folder selection
        output_layout = QGridLayout()
        output_layout.addWidget(QLabel('Select Output Folder:'), 0, 0)
        self.output_path_label = QLabel()
        output_layout.addWidget(self.output_path_label, 1, 0, 1, 2)
        output_select_button = QPushButton('Select', self)
        output_select_button.clicked.connect(self.select_output_folder)
        output_layout.addWidget(output_select_button, 0, 1)
        layout.addLayout(output_layout)

        # Run and Cancel buttons
        buttons_layout = QHBoxLayout()
        run_button = QPushButton('Run', self)
        run_button.clicked.connect(self.run_analysis)
        buttons_layout.addWidget(run_button)
        cancel_button = QPushButton('Cancel', self)
        cancel_button.clicked.connect(self.cancel_analysis)
        buttons_layout.addWidget(cancel_button)
        layout.addLayout(buttons_layout)

        self.setLayout(layout)

        # Connect signals
        self.signal.started.connect(self.show_running_msg)
        self.signal.completed.connect(self.show_complete_msg)
        self.signal.error.connect(self.show_error_msg)

    def select_forward_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, 'Select Forward Files', '', 'FASTQ Files (*.fastq *.fq *fastq.gz *fq.gz)'
        )
        self.forward_files.extend(file_paths)
        self.update_forward_paths()

    def add_forward_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, 'Add Forward Files', '', 'FASTQ Files (*.fastq *.fq *fastq.gz *fq.gz)'
        )
        self.forward_files.extend(file_paths)
        self.update_forward_paths()

    def delete_forward_file(self):
        selected_index = self.forward_paths.textCursor().blockNumber()
        if selected_index >= 0 and selected_index < len(self.forward_files):
            del self.forward_files[selected_index]
            self.update_forward_paths()

    def clear_forward_files(self):
        self.forward_files.clear()
        self.update_forward_paths()

    def select_reverse_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, 'Select Reverse Files', '', 'FASTQ Files (*.fastq *.fq *fastq.gz *fq.gz)'
        )
        self.reverse_files.extend(file_paths)
        self.update_reverse_paths()

    def add_reverse_files(self):
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, 'Add Reverse Files', '', 'FASTQ Files (*.fastq *.fq *fastq.gz *fq.gz)'
        )
        self.reverse_files.extend(file_paths)
        self.update_reverse_paths()

    def delete_reverse_file(self):
        selected_index = self.reverse_paths.textCursor().blockNumber()
        if selected_index >= 0 and selected_index < len(self.reverse_files):
            del self.reverse_files[selected_index]
            self.update_reverse_paths()

    def clear_reverse_files(self):
        self.reverse_files.clear()
        self.update_reverse_paths()

    def select_output_folder(self):
        output_folder = QFileDialog.getExistingDirectory(self, 'Select Output Folder')
        if output_folder:
            self.output_folder = output_folder
            self.update_output_path()

    def update_forward_paths(self):
        paths = '\n'.join(self.forward_files)
        self.forward_paths.setPlainText(paths)

    def update_reverse_paths(self):
        paths = '\n'.join(self.reverse_files)
        self.reverse_paths.setPlainText(paths)

    def update_output_path(self):
        self.output_path_label.setText(self.output_folder)

    def run_analysis(self):
        if not self.forward_files or not self.reverse_files or not self.output_folder:
            QMessageBox.critical(
                self, 'Error', 'Please select both forward and reverse files, and an output folder.'
            )
            return

        # Command to be executed within the conda environment
        command = ['Trinity', '--seqType', 'fq']
        command.extend(['--left', ','.join(self.forward_files)])
        command.extend(['--right', ','.join(self.reverse_files)])
        command.extend([
            '--max_memory', '50G',
                        '--output', f'{self.output_folder}/trinity_assembly',
            '--CPU', '16'
        ])

        # Start the timer
        start_time = datetime.now()

        def run_trinity():
            self.signal.started.emit()

            try:
                # Activate the trinity_env environment and run the command
                activate_command = ['conda', 'run', '-n', 'trinity_env'] + command
                process = subprocess.Popen(activate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()

                if process.returncode == 0:
                    # Move the Trinity assembly files to the destination folder
                    source_folder = self.output_folder
                    destination_folder = os.path.join(self.output_folder, "trinity_assembly")
                    file_extensions = [".Trinity.fasta", ".Trinity.fasta.gene_trans_map"]

                    for file_name in os.listdir(source_folder):
                        if any(file_name.endswith(extension) for extension in file_extensions):
                            source_path = os.path.join(source_folder, file_name)
                            destination_path = os.path.join(destination_folder, file_name)
                            shutil.move(source_path, destination_path)
                    end_time = datetime.now()  # End time
                    elapsed_time = end_time - start_time  # Elapsed time
                    elapsed_time_str = str(elapsed_time)

                    if process.returncode == 0:
                        self.signal.completed.emit(elapsed_time.total_seconds())
                    else:
                        self.signal.error.emit("Trinity command failed. Please check your inputs and try again.")
                else:
                    self.signal.error.emit("Trinity command failed. Please check your inputs and try again.")
            except Exception as e:
                self.signal.error.emit(f"An error occurred while running Trinity: {str(e)}")

        # Start the thread
        threading.Thread(target=run_trinity, daemon=True).start()

    @pyqtSlot()
    def show_running_msg(self):
        self.running_msg = QMessageBox(self)
        self.running_msg.setIcon(QMessageBox.Information)
        self.running_msg.setText(" ")
        self.running_msg.setWindowTitle("Running")
        self.running_msg.setStandardButtons(QMessageBox.Cancel)
        self.running_msg.buttonClicked.connect(self.cancel_analysis)
        self.running_msg.show()

    @pyqtSlot(float)
    def show_complete_msg(self, duration):
        if self.running_msg:
            self.running_msg.close()  # Close the "Running" message box

        # Convert duration to hours, minutes, and seconds
        hours = int(duration / 3600)
        minutes = int((duration % 3600) / 60)
        seconds = int(duration % 60)

        # Format the duration as "hours minutes seconds"
        duration_str = f"{hours}h {minutes}m {seconds}s"

        # Show "Assembly is complete" message box with timestamp
        complete_msg = QMessageBox(self)
        complete_msg.setIcon(QMessageBox.Information)
        complete_msg.setText(f"Assembly is completed.\nDuration: {duration_str}")
        complete_msg.setWindowTitle("Complete")
        complete_msg.setStandardButtons(QMessageBox.Ok)
        complete_msg.buttonClicked.connect(self.close_windows)
        complete_msg.show()

    @pyqtSlot(str)
    def show_error_msg(self, error_message):
        if self.running_msg:
            self.running_msg.close()  # Close the "Running" message box

        QMessageBox.critical(self, "Trinity Error", error_message)

    def cancel_analysis(self):
        if self.process and self.process.poll() is None:
            # Process is still running, terminate it
            self.process.terminate()
            self.process.wait()

        if self.running_msg:
            self.running_msg.close()  # Close the "Running" message box

    def close_windows(self):
        if self.running_msg:
            self.running_msg.close()  # Close the "Running" message box
        self.close()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = PopupWindow()
    window.show()
    sys.exit(app.exec_())