import grpc
from concurrent import futures
import time
import example_pb2
import example_pb2_grpc

# 实现Greeter服务
class GreeterServicer(example_pb2_grpc.GreeterServicer):
    def SayHello(self, request, context):
        return example_pb2.Response(message=f"Hello, {request.name}!")

# 启动gRPC服务器
def serve():
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=96))
    example_pb2_grpc.add_GreeterServicer_to_server(GreeterServicer(), server)
    server.add_insecure_port('[::]:50051')
    print("Server started on port 50051")
    server.start()
    server.wait_for_termination()

if __name__ == '__main__':
    serve()
