import grpc
import time
import example_pb2
import example_pb2_grpc
import multiprocessing

# 定义每个客户端进程发送请求并测量延迟
def send_requests(client, name, num_requests):
    total_delay = 0.0
    
    # 每个客户端进程发送多个请求
    for _ in range(num_requests):
        start_time = time.time()
        response = client.SayHello(example_pb2.Request(name=name))
        # local_hello()
        end_time = time.time()
        
        # 计算延迟并累计
        delay = end_time - start_time
        total_delay += delay
        # print(f"Request to {name} took {delay:.4f} seconds.")
    
    # 返回当前进程的平均延迟
    avg_delay = total_delay / num_requests
    return avg_delay

# 用于创建gRPC通道和客户端
def create_client():
    channel = grpc.insecure_channel('localhost:50051')
    return example_pb2_grpc.GreeterStub(channel)

def local_hello():
    a,b= 1,2
    c= a+b

# 创建多个进程，模拟多个客户端并发请求
def run(num_processes, num_requests_per_process):
    processes = []
    results = []

    start_time = time.time()
    # 创建多个进程，每个进程发送多个请求
    for i in range(num_processes):
        process = multiprocessing.Process(target=worker, args=(i, num_requests_per_process, results))
        processes.append(process)
        process.start()

    # 等待所有进程完成
    for process in processes:
        process.join()

    end_time= time.time()
    print(f"Overall Throughput: {num_processes * num_requests_per_process / (end_time - start_time)} requests per second")
    # 计算所有进程的平均延迟
    # avg_delay_all = sum(results) / len(results)
    # print(f"Average delay for all requests: {avg_delay_all:.4f} seconds.")

# 每个工作进程的任务
def worker(process_id, num_requests, results):
    client = create_client()
    avg_delay = send_requests(client, f"Client-{process_id}", num_requests)
    results.append(avg_delay)
    print(f"Average delay for Client-{process_id}: {avg_delay*1000000:.4f} us.")

if __name__ == '__main__':
    num_processes = 3
    num_requests_per_process = 1000  # 每个进程发送10个请求
    run(num_processes, num_requests_per_process)
