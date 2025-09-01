#!/usr/bin/env python3
"""
Test script for admin endpoints
"""
import requests
import json
from datetime import datetime

BASE_URL = "http://localhost:8000/api"

def test_admin_endpoints():
    print("Testing Admin Endpoints\n")
    print("=" * 50)
    
    # Step 1: Login as admin to get access token
    print("\n1. Testing Admin Login...")
    login_data = {
        "email": "admin@example.com",  # Replace with your admin credentials
        "password": "adminpassword"
    }
    
    try:
        response = requests.post(f"{BASE_URL}/auth/login/", json=login_data)
        
        if response.status_code == 200:
            data = response.json()
            access_token = data.get("accessToken")
            user_info = data.get("user")
            
            if not user_info.get("is_staff"):
                print("✗ User is not an admin/staff member")
                return
                
            print(f"✓ Admin login successful!")
            print(f"  User: {user_info.get('email')}")
            print(f"  Is Staff: {user_info.get('is_staff')}")
            print(f"  Is Superuser: {user_info.get('is_superuser')}")
            
            headers = {"Authorization": f"Bearer {access_token}"}
            
            # Step 2: Test Dashboard Stats
            print("\n2. Testing Dashboard Stats...")
            stats_response = requests.get(f"{BASE_URL}/admin/dashboard/stats/", headers=headers)
            
            if stats_response.status_code == 200:
                stats = stats_response.json()
                print(f"✓ Dashboard stats retrieved successfully!")
                print(f"  Total Users: {stats['user_stats']['total']}")
                print(f"  Active Users: {stats['user_stats']['active']}")
                print(f"  Total Conversations: {stats['conversation_stats']['total']}")
                print(f"  Countries Reached: {len(stats['users_by_country'])}")
            else:
                print(f"✗ Failed to get dashboard stats: {stats_response.status_code}")
                print(f"  Response: {stats_response.text}")
            
            # Step 3: Test Get Users
            print("\n3. Testing Get Users...")
            users_response = requests.get(
                f"{BASE_URL}/admin/users/",
                headers=headers,
                params={"page": 1, "page_size": 10}
            )
            
            if users_response.status_code == 200:
                users_data = users_response.json()
                print(f"✓ Users retrieved successfully!")
                print(f"  Total Users: {users_data['total']}")
                print(f"  Page: {users_data['page']}/{users_data['total_pages']}")
                print(f"  Users on this page: {len(users_data['users'])}")
                
                if users_data['users']:
                    # Test user detail endpoint with first user
                    test_user_id = users_data['users'][0]['id']
                    
                    # Step 4: Test Get User Detail
                    print(f"\n4. Testing Get User Detail (ID: {test_user_id})...")
                    detail_response = requests.get(
                        f"{BASE_URL}/admin/users/{test_user_id}/",
                        headers=headers
                    )
                    
                    if detail_response.status_code == 200:
                        user_detail = detail_response.json()
                        print(f"✓ User detail retrieved successfully!")
                        print(f"  Username: {user_detail['username']}")
                        print(f"  Email: {user_detail['email']}")
                        print(f"  Total Conversations: {user_detail.get('total_conversations', 0)}")
                    else:
                        print(f"✗ Failed to get user detail: {detail_response.status_code}")
            else:
                print(f"✗ Failed to get users: {users_response.status_code}")
                print(f"  Response: {users_response.text}")
            
            # Step 5: Test Search Users
            print("\n5. Testing User Search...")
            search_response = requests.get(
                f"{BASE_URL}/admin/users/",
                headers=headers,
                params={"search": "test", "page": 1, "page_size": 10}
            )
            
            if search_response.status_code == 200:
                search_data = search_response.json()
                print(f"✓ User search successful!")
                print(f"  Found {search_data['total']} users matching 'test'")
            else:
                print(f"✗ Failed to search users: {search_response.status_code}")
            
            # Step 6: Test Export Users CSV
            print("\n6. Testing Export Users CSV...")
            export_response = requests.get(f"{BASE_URL}/admin/users/export/", headers=headers)
            
            if export_response.status_code == 200:
                print(f"✓ CSV export successful!")
                print(f"  Response size: {len(export_response.content)} bytes")
                print(f"  Content type: {export_response.headers.get('content-type')}")
                
                # Save CSV to file for inspection
                filename = f"users_export_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
                with open(filename, 'wb') as f:
                    f.write(export_response.content)
                print(f"  Saved to: {filename}")
            else:
                print(f"✗ Failed to export CSV: {export_response.status_code}")
                
        else:
            print(f"✗ Login failed: {response.status_code}")
            print(f"  Response: {response.text}")
            print("\nNote: You need to create an admin user first or use valid admin credentials")
            print("To create an admin user, run:")
            print("  python manage.py createsuperuser")
            
    except requests.exceptions.ConnectionError:
        print("✗ Could not connect to server. Make sure the Django server is running on port 8000")
    except Exception as e:
        print(f"✗ Error: {e}")

if __name__ == "__main__":
    test_admin_endpoints()
    print("\n" + "=" * 50)
    print("Test completed!")
