./a.out &  # Replace "my_program" with the actual program name or path
PROGRAM_PID=$!

# Check if the program started successfully
if [ -z "$PROGRAM_PID" ]; then
  echo "Failed to start the program."
  exit 1
fi

echo "Program started with PID: $PROGRAM_PID"

# Start monitoring CPU usage with pidstat for the specific PID
# Adjust the interval (1 second here) and output file name if needed
echo "Starting CPU monitoring for PID $PROGRAM_PID"
pidstat -p $PROGRAM_PID 1 > cpu_usage_$PROGRAM_PID.log &

# Wait for the program to finish
wait $PROGRAM_PID

# Stop monitoring once the program completes
echo "Program completed. CPU monitoring stopped."
